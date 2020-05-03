
export tfullcontract, fullcontract

function fullcontract(contractor, region::T, tol = default_tolerance, maxiters = 100) where T
    for i in 1:maxiters
        region2 = region .∩ contractor(region)
        if maximum(diam.(region2)) < maximum(diam.(region))
            region = region2
        else
            break
        end
        maximum(diam.(region)) < tol && break
    end
    return region
end

tfullcontract(contractor, regions, tol = default_tolerance) =
    qmap(r -> fullcontract(contractor, r, tol), collect(regions))
