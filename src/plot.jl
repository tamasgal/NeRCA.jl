"""
    function savefigs(basename::AbstractString; types=["svg", "tex", "pdf", "png"])

Save the current plot in different formats.
"""
function savefigs(basename::AbstractString;
                  types=["svg", "tex", "pdf", "png"])
    print("Saving ")
    prefix = ""
    suffix = ", "
    for (idx, filetype) in enumerate(types)
        filename = basename * '.' * filetype
        if idx == length(types) - 1
            suffix = ""
        end
        if idx == length(types)
            prefix = " and "
            suffix = "."
        end
        print(prefix * filename * suffix)
        savefig(filename)
    end
end
