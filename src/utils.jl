function float_string(x, n, sym = "p")
    order = Int(floor(log10(x) + 1))
    str = rpad(round(x, digits=n), 4 + order, "0")
    str[1:order] * sym * str[order + 2:end]
end

function ExportToMathematicaInterp(array::Array{Float64,2}, name)
    ret = string(name) * " = {"

    for i = 1:size(array)[1]
        for j = 1:size(array)[2]
            ret = ret * "{" * string(x[i, 1, 1]) * ", " *
                               string(z[1, 1, j]) * ", " *
                       replace(string(array[i, j]), "e" => "*^") * "}"
	     if i*j < size(array)[1]*size(array)[2]
	        ret = ret * ", "
	    end
        end
    end

    ret = ret * "}\n\n"

    ret
end
