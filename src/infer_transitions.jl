function infer_transitions(possible_transition)
  nstates=size(possible_transition,1)
  
  possible_transition_tuple = Vector{Tuple{Int64,Int64}}(undef,0)
  for idx_begin = 1:nstates
    for idx_end = 1:nstates
      if possible_transition[idx_begin,idx_end]
        push!(possible_transition_tuple,CartesianIndex(idx_begin,idx_end))
      end
    end
  end

  reachable_within=fill(Inf,nstates,nstates)
  for ii in nstates:(-1):0
      for idx in findall((possible_transition^ii).>=1)
          reachable_within[idx]=ii
      end
  end
  possible_sub=.!(isinf.(reachable_within))
  
  F_next_states=Vector{Vector{Int64}}(undef, nstates)
  F_sub_states=Matrix{Vector{Int64}}(undef, nstates, nstates)
  for idx_begin in 1:nstates
    idx_sub=possible_transition[idx_begin,:]
    if any(idx_sub) F_next_states[idx_begin]=findall(idx_sub) end
  end
  for idx_begin in 1:nstates
    for idx_end in 1:nstates
      idx_sub=possible_sub[idx_begin,:].&possible_sub[:,idx_end]
      if any(idx_sub) F_sub_states[idx_begin,idx_end]=findall(idx_sub) end
    end
  end

  absorbing_states=.!any(possible_transition, dims=2)
  
  cycles=fill(false,nstates)
  for ii in nstates:(-1):1
    for idx in findall(LinearAlgebra.diag(possible_transition^ii).>=1)
      cycles[idx]=true
    end
  end

  return (
    possible_transition_tuple,
    reachable_within,
    possible_sub,
    F_next_states,
    F_sub_states,
    absorbing_states,
    cycles)
end