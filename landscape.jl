using Graphs


# Create a graph to represent locations where
# directed weights indicate the hazard rate for infection
# from one location to the next.

# rectangular bounds [ left, right, up, down]
function poisson_point_process(rectangular_bounds, cnt)
    left, right, down, up=rectangular_bounds
    width_scale=right-left
    length_scale=up-down
    x=zeros(Float64, 2, cnt)
    for draw_idx in 1:cnt
    	x[:,draw_idx]=[left+rand()*width_scale, down+rand()*length_scale]
    end
    x
end


# rectangular bounds [ left, right, up, down]
# pi x (hard sphere radius^2)/area=areal_fraction
# Safe to keep this below 0.5. They get packed in there.
function hard_sphere_process(rectangular_bounds, cnt, areal_fraction)
    left, right, down, up=rectangular_bounds
    width_scale=right-left
    length_scale=up-down

    min_distance2=4*areal_fraction*width_scale*length_scale/(pi*cnt)

    x=zeros(Float64, 2, cnt)
    draw_idx=1
    found_idx=0
    while found_idx<cnt && draw_idx<cnt*100000
        xy=[left+rand()*width_scale, down+rand()*length_scale]
        found=true
        for search_idx in 1:found_idx
            distance2=(x[1,search_idx]-xy[1])^2 + (x[2,search_idx]-xy[2])^2
            if distance2<min_distance2
                found=false
                break
            end
        end
        if found
            found_idx+=1
            x[:,found_idx]=xy
        end
        draw_idx+=1
    end
    x
end



function all_to_all_adjacency_matrix(locations)
    cnt=size(locations)[2]
    adjacency_matrix=zeros(Float64, cnt, cnt)
    for i in 1:cnt
    	for j in (i+1):cnt
    		distance_squared=((locations[1,i]-locations[1,j])^2 +
    			(locations[2,i]-locations[2,j])^2)
    		adjacency_matrix[i,j]=distance_squared
    		adjacency_matrix[j,i]=distance_squared
    	end
    end
    adjacency_matrix
end


# Given all neighbors of a location, assign propensities for
# transmission to that neighbor based on distance.
# This gives a larger total propensity if more neighbors are closer.
function hop_jump_kernel(distance_matrix, cutoff, hop_propensity, jump_propensity)
    propensity_matrix=zeros(Float64, size(distance_matrix)...)
    cnt=size(distance_matrix)[1]
    for from_idx in 1:cnt
        for to_idx in 1:cnt
        	if from_idx!=to_idx
	        	if (distance_matrix[from_idx, to_idx]<cutoff)
	        	    propensity_matrix[from_idx, to_idx]=hop_propensity
	        	else
	        	    propensity_matrix[from_idx, to_idx]=jump_propensity
	        	end
	        end
        end
    end
    propensity_matrix
end



function hop_jump_kernel_normalized(distance_matrix, cutoff, base_propensity,
       total_propensity)
	propensity_matrix=hop_jump_kernel(distance_matrix, cutoff, base_propensity)
    cnt=size(distance_matrix)[1]
	for from_idx in 1:cnt
	    normalization=sum(propensity_matrix[from_idx, :])
        propensity_matrix[from_idx,:]*=total_propensity/normalization
	end
    propensity_matrix
end



