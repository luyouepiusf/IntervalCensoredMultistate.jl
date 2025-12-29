function reg_est_bone(
  TTij,
  ssij,
  nobss,
  zzi,
  F_zidx,
  possible_transition,
  ntimepoints;
  niter=1000,tol=1e-4,factor=0.8)
  
  nind=size(TTij,1)
  ndimz=size(zzi,2)
  nstates=size(possible_transition,1)
  eye_nstates=Matrix{Float64}(LinearAlgebra.I,nstates,nstates)
  (
    possible_transition_tuple,
    reachable_within,
    possible_sub,
    F_next_states,
    F_sub_states,
    absorbing_states,
    cycles)=infer_transitions(possible_transition)
  
  F_nzidx=fill(0,nstates,nstates)
  for (idx_begin,idx_end) in possible_transition_tuple
    F_nzidx[idx_begin,idx_end]=length(F_zidx[idx_begin,idx_end])
  end
  
  F_betaz_est = Matrix{Vector}(undef, nstates, nstates)
  F_betaz_old = Matrix{Vector}(undef, nstates, nstates)
  F_eps = fill(0.0, nstates, nstates)
  F_hazard0 = Matrix{Vector}(undef, nstates, nstates)
  
  for idx_begin = 1:nstates
    for idx_end = 1:nstates
      if possible_transition[idx_begin,idx_end]
        F_betaz_est[idx_begin,idx_end] = fill(0.0, F_nzidx[idx_begin,idx_end])
        F_betaz_old[idx_begin,idx_end] = fill(0.0, F_nzidx[idx_begin,idx_end])
        F_hazard0[idx_begin,idx_end] = fill(0.002, ntimepoints)
      end
    end
  end
  F_SS_betaz = Matrix{Vector}(undef, nstates, nstates)
  F_II_betaz = Matrix{Matrix}(undef, nstates, nstates)

  F_nonzero_g=Matrix{Vector{Int64}}(undef,nind,size(ssij,2))
  F_nonzero_gd=Array{Vector{Int64},3}(undef,nind,size(ssij,2),nstates)
  for ii in 1:nind
    for jj in 1:(nobss[ii]-1)
      a_sub_states=Vector{Int64}(undef,0)
      for idx_begin in findall(ssij[ii,jj,:])
        for idx_end in findall(ssij[ii,jj+1,:])
          if isassigned(F_sub_states,idx_begin,idx_end) append!(a_sub_states,F_sub_states[idx_begin,idx_end]) end
        end
      end
      unique_sub_states=sort(unique(a_sub_states))
      F_nonzero_g[ii,jj]=unique_sub_states
      a_next_states=Vector{Int64}(undef,0)
      for idx_begin in unique_sub_states
        if isassigned(F_next_states,idx_begin) 
          F_nonzero_gd[ii,jj,idx_begin]=F_next_states[idx_begin] 
        else
          F_nonzero_gd[ii,jj,idx_begin]=[]
        end
      end
    end
  end

  F_one_g=fill(false,nind,size(ssij,2),nstates)
  for ii in 1:nind
    for jj in 1:(nobss[ii]-1)
      if sum(ssij[ii,jj,:])==1 && all(ssij[ii,jj,:].==ssij[ii,jj+1,:])
        ss=findfirst(ssij[ii,jj,:])
        if !cycles[ss]
          F_one_g[ii,jj,ss]=true
        end
      end
    end
  end

  ## define blocks
  nblock=fill(0,nind)
  block2jj=fill(-1,nind,size(TTij)[2],2)
  block2tt=fill(-1,nind,size(TTij)[2],2)
  for ii in 1:nind
    temp_nblock=0
    temp_block2jj=fill(-1,size(TTij)[2],2)
    temp_block2tt=fill(-1,size(TTij)[2],2)
    temp_begin=nobss[ii]
    temp_end=nobss[ii]
    for jj in (nobss[ii]-1):(-1):1
      if sum(ssij[ii,jj,:])==1
        temp_begin=jj
        temp_nblock=temp_nblock+1
        temp_block2jj[temp_nblock,1]=temp_begin
        temp_block2jj[temp_nblock,2]=temp_end
        temp_block2tt[temp_nblock,1]=TTij[ii,temp_begin]
        temp_block2tt[temp_nblock,2]=TTij[ii,temp_end]
        temp_end=jj
      end
    end
    nblock[ii]=temp_nblock
    block2jj[ii,1:temp_nblock,:]=temp_block2jj[temp_nblock:(-1):1,:]
    block2tt[ii,1:temp_nblock,:]=temp_block2tt[temp_nblock:(-1):1,:]
  end

  prob_by_jj=fill(NaN,nstates,nstates,size(TTij)[2])
  prob_lt_tt=fill(NaN,nstates,nstates,ntimepoints)
  prob_gt_tt=fill(NaN,nstates,nstates,ntimepoints)
  prob_ge_tt=fill(NaN,nstates,nstates,ntimepoints)

  F_Egd = fill(0.0, nind, nstates, nstates, ntimepoints)
  F_Eg = fill(0.0, nind, nstates, ntimepoints)
  F_expr = fill(NaN, nind, nstates, nstates)
  F_zexpr = Array{Vector{Float64},3}(undef, nind, nstates, nstates)
  F_zzexpr = Array{Matrix{Float64},3}(undef, nind, nstates, nstates)

  start_time=time()
  for iter in 1:niter
    println(iter," - ")
    current_time=time()
    println(current_time-start_time)
    for ii in 1:nind
      prob=zeros(nstates,nstates,ntimepoints)
      expr=zeros(nstates,nstates)
      for (idx_begin,idx_end) in possible_transition_tuple
        zz_temp=zzi[ii,F_zidx[idx_begin,idx_end]]
        expr[idx_begin,idx_end]=exp(LinearAlgebra.dot(F_betaz_est[idx_begin,idx_end],zz_temp))
        F_expr[ii,idx_begin,idx_end]=expr[idx_begin,idx_end]
        F_zexpr[ii,idx_begin,idx_end] = F_expr[ii,idx_begin,idx_end] * zz_temp
        F_zzexpr[ii,idx_begin,idx_end] = F_zexpr[ii,idx_begin,idx_end] * transpose(zz_temp)
        prob[idx_begin,idx_end,:]=1.0.-exp.(-F_hazard0[idx_begin,idx_end]*expr[idx_begin,idx_end])
      end
      for idx_begin in 1:nstates
        prob[idx_begin,idx_begin,:]=1.0.-sum(prob[idx_begin,:,:],dims=1)
      end
      prob_by_block=fill(NaN,size(ssij)[2]) 
      prob_by_jj=fill(NaN,nstates,nstates,size(ssij)[2])
      prob_lt_tt=fill(NaN,nstates,nstates,ntimepoints) # can be defined outside the loop
      prob_gt_tt=fill(NaN,nstates,nstates,ntimepoints) # can be defined outside the loop
      prob_ge_tt=fill(NaN,nstates,nstates,ntimepoints) # can be defined outside the loop

      for jj in 1:(nobss[ii]-1)
        prob_temp=copy(eye_nstates)
        for tt in (TTij[ii,jj]+1):TTij[ii,jj+1]
          prob_lt_tt[:,:,tt]=prob_temp
          prob_temp=prob_temp*prob[:,:,tt]
        end
        prob_by_jj[:,:,jj]=prob_temp
        prob_temp=copy(eye_nstates)
        for tt in (TTij[ii,jj+1]):(-1):(TTij[ii,jj]+1)
          prob_gt_tt[:,:,tt]=prob_temp
          prob_temp=prob[:,:,tt]*prob_temp
          prob_ge_tt[:,:,tt]=prob_temp
        end
      end
      F_lt_jj = Vector{Matrix{Float64}}(undef,size(ssij)[2])
      F_gt_jj = Vector{Matrix{Float64}}(undef,size(ssij)[2])
      for bb in 1:nblock[ii]
        block_temp=fill(1.0,1,sum(ssij[ii,block2jj[ii,bb,1],:]))
        for jj in block2jj[ii,bb,1]:(block2jj[ii,bb,2]-1)
          F_lt_jj[jj]=block_temp
          block_temp=block_temp*prob_by_jj[ssij[ii,jj,:],ssij[ii,jj+1,:],jj]
        end
        prob_by_block[bb]=sum(block_temp)
        block_temp=fill(1.0,sum(ssij[ii,block2jj[ii,bb,2],:]),1)
        for jj in (block2jj[ii,bb,2]-1):(-1):block2jj[ii,bb,1]
          F_gt_jj[jj]=block_temp
          block_temp=prob_by_jj[ssij[ii,jj,:],ssij[ii,jj+1,:],jj]*block_temp
        end
      end
      for bb in 1:nblock[ii]
        for jj in block2jj[ii,bb,1]:(block2jj[ii,bb,2]-1)
          for tt in (TTij[ii,jj]+1):(TTij[ii,jj+1])
            for a_sub_idx_begin in F_nonzero_g[ii,jj]
              if F_one_g[ii,jj,a_sub_idx_begin]
                g_temp=1.0
              else
                g_temp=sum(
                  F_lt_jj[jj]*
                  prob_lt_tt[ssij[ii,jj,:],a_sub_idx_begin:a_sub_idx_begin,tt]*
                  prob_ge_tt[a_sub_idx_begin:a_sub_idx_begin,ssij[ii,jj+1,:],tt]*
                  F_gt_jj[jj])/prob_by_block[bb]
              end
              F_Eg[ii,a_sub_idx_begin,tt]=g_temp
              for a_sub_idx_end in F_nonzero_gd[ii,jj,a_sub_idx_begin]
                gd_temp=sum(
                  F_lt_jj[jj]*
                  prob_lt_tt[ssij[ii,jj,:],a_sub_idx_begin:a_sub_idx_begin,tt]*
                  prob[a_sub_idx_begin:a_sub_idx_begin,a_sub_idx_end:a_sub_idx_end,tt]*
                  prob_gt_tt[a_sub_idx_end:a_sub_idx_end,ssij[ii,jj+1,:],tt]*
                  F_gt_jj[jj])/prob_by_block[bb]
                F_Egd[ii,a_sub_idx_begin,a_sub_idx_end,tt]=gd_temp
              end
            end
          end
        end
      end
    end
    for (idx_begin,idx_end) in possible_transition_tuple
      AA0=fill(0.0,ntimepoints)
      AA1=fill(0.0,F_nzidx[idx_begin,idx_end],ntimepoints)
      AA2=fill(0.0,F_nzidx[idx_begin,idx_end],F_nzidx[idx_begin,idx_end],ntimepoints)
      BB1=fill(0.0,F_nzidx[idx_begin,idx_end],ntimepoints)
      BB2=fill(0.0,F_nzidx[idx_begin,idx_end],F_nzidx[idx_begin,idx_end],ntimepoints)
      for tt in 1:ntimepoints
        for ii in 1:nind
          AA0[tt]=AA0[tt]+(F_Eg[ii,idx_begin,tt]-F_Egd[ii,idx_begin,idx_end,tt]/2.0)*F_expr[ii,idx_begin,idx_end]
          AA1[:,tt]=AA1[:,tt]+(F_Eg[ii,idx_begin,tt]-F_Egd[ii,idx_begin,idx_end,tt]/2.0)*F_zexpr[ii,idx_begin,idx_end]
          AA2[:,:,tt]=AA2[:,:,tt]+(F_Eg[ii,idx_begin,tt]-F_Egd[ii,idx_begin,idx_end,tt]/2.0)*F_zzexpr[ii,idx_begin,idx_end]
        end
        if AA0[tt]==0.0
          BB1[:,tt].=0.0
          BB2[:,:,tt].=0.0
        else
          BB1[:,tt]=AA1[:,tt]/AA0[tt]
          BB2[:,:,tt]=AA2[:,:,tt]/AA0[tt]
        end
      end
      SS_betaz=fill(0.0,F_nzidx[idx_begin,idx_end]) # can be defined outside the loop
      II_betaz=fill(0.0,F_nzidx[idx_begin,idx_end],F_nzidx[idx_begin,idx_end]) # can be defined outside the loop
      for tt in 1:ntimepoints
        for ii in 1:nind
          SS_betaz=SS_betaz+F_Egd[ii,idx_begin,idx_end,tt]*zzi[ii,F_zidx[idx_begin,idx_end]]-F_Egd[ii,idx_begin,idx_end,tt]*BB1[:,tt]
          II_betaz=II_betaz+F_Egd[ii,idx_begin,idx_end,tt]*(BB1[:,tt]*transpose(BB1[:,tt])-BB2[:,:,tt])
        end
      end
      F_betaz_old[idx_begin,idx_end]=deepcopy(F_betaz_est[idx_begin,idx_end])
      F_betaz_est[idx_begin,idx_end]=F_betaz_est[idx_begin,idx_end]-factor*(II_betaz\SS_betaz)
      F_eps[idx_begin,idx_end]=sum(abs.(F_betaz_est[idx_begin,idx_end]-F_betaz_old[idx_begin,idx_end]))
      F_II_betaz[idx_begin,idx_end]=II_betaz
      for tt = 1:ntimepoints
        sum_Egd = sum(F_Egd[:,idx_begin,idx_end,tt]) # can be defined outside the loop
        if AA0[tt]==0.0
          F_hazard0[idx_begin,idx_end][tt]=0.0
        else
          F_hazard0[idx_begin,idx_end][tt]=sum_Egd/AA0[tt]
        end
      end
    end
    max_eps=sum(F_eps)
    if max_eps<tol break end
  end
  
  return(
    F_betaz_est,
    F_II_betaz,
    F_hazard0,
    F_Eg,
    F_eps)
end

