function simplify(
  ssij,TTij,nobss,
  possible_transition;
  less_time=true,
  non_info=false)

  (
    possible_transition_tuple,
    reachable_within,
    possible_sub,
    F_next_states,
    F_sub_states,
    absorbing_states,
    cycles)=infer_transitions(possible_transition)

  ssij1=fill(false,size(ssij))
  TTij1=fill(-100,size(TTij))
  nobss1=fill(-100,size(nobss))

  for ii in 1:size(TTij,1)
    temp_nobss=nobss[ii]
    temp_TTij=TTij[ii,1:temp_nobss]
    temp_ssij=ssij[ii,1:temp_nobss,:]
    for iter in 1:3
      for jj in 2:temp_nobss
        temp_from=temp_ssij[jj-1,:]
        temp_to=temp_ssij[jj,:]
        for idx_begin in findall(temp_from)
          if all(possible_sub[idx_begin,temp_to].==false)
            temp_ssij[jj-1,idx_begin]=false
          end
        end
        for idx_end in findall(temp_to)
          if all(possible_sub[temp_from,idx_end].==false)
            temp_ssij[jj,idx_end]=false
          end
        end
      end
    end
    if less_time
      # single state
      idx=fill(true,length(temp_TTij))
      for jj in 1:length(temp_TTij)
        if jj==1 continue end
        if jj==length(temp_TTij) continue end
        if sum(temp_ssij[jj-1,:])==1 && sum(temp_ssij[jj,:])==1 && sum(temp_ssij[jj+1,:])==1
          ss_jjm1=findfirst(temp_ssij[jj-1,:])
          ss_jj=findfirst(temp_ssij[jj,:])
          ss_jjp1=findfirst(temp_ssij[jj+1,:])
          if ss_jjm1==ss_jj && ss_jj==ss_jjp1 && !cycles[ss_jj]
            idx[jj]=false 
          end
        end
      end
      temp_ssij=temp_ssij[idx,:]
      temp_TTij=temp_TTij[idx]
      
      # end state
      idx=fill(true,length(temp_TTij))
      for jj in 1:length(temp_TTij)
        if jj==1 continue end
        if sum(temp_ssij[jj,:])==1
          if all(absorbing_states[temp_ssij[jj,:]]) 
            idx[(jj+1):end].=false
            break
          end
        end
      end
      temp_ssij=temp_ssij[idx,:]
      temp_TTij=temp_TTij[idx]

      # noninformative state
      if non_info
        idx=fill(true,length(temp_TTij))
        for jj in 1:length(temp_TTij)
          if jj==1 continue end
          temp_possible_future=map(any,eachcol(possible_transition[temp_ssij[jj-1,:],:]))
          if all(temp_ssij[jj,temp_possible_future])
            idx[jj]=false
          end
        end
        temp_ssij=temp_ssij[idx,:]
        temp_TTij=temp_TTij[idx]
      end
    end
    TTij1[ii,1:length(temp_TTij)]=temp_TTij
    ssij1[ii,1:length(temp_TTij),:]=temp_ssij
    nobss1[ii]=length(temp_TTij)
  end
  nmaxs=maximum(nobss1)
  return(ssij1[:,1:nmaxs,:],TTij1[:,1:nmaxs],nobss1)
end
