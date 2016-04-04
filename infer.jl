
function single_trajectory(cnt)
    params=Dict([("gamma", 0.15), ("betahop", 0.5), ("betajump", 0.05), ("cutoff", 0.3)])
    points=poisson_point_process([0, 1, 0, 1], cnt)
    distances=all_to_all_adjacency_matrix(points)
    propensity_matrix=hop_jump_kernel(distances, params["cutoff"], params["betahop"],
        params["betajump"])

    rand_seed=9324
    sir_times=single_hop_skip(params, propensity_matrix, rand_seed)
end


sir_times=single_trajectory(50)

println(size(sir_times))
println(sir_times)

