src_path = raw"C:\Users\smmywael\.julia\dev\OpPoDyn\test\WECC_model_tests\PV_pf\variables-testcase3Bus-with-event_SystemBase_IdealSlack.csv"
main_path = raw"C:\Users\smmywael\.julia\dev\OpPoDyn\test\WECC_model_tests\PV_pf\variables-testcase3Bus-with-event_SystemBase_IdealSlack.csv"

function parse_row(line)
    cells = split(rstrip(line, ['\n', '\r']), ';')
    return [strip(c, '"') for c in cells]
end

function read_first_rows(path, n)
    rows = Vector{Vector{String}}()
    open(path, "r") do f
        for _ in 1:n
            line = readline(f)
            push!(rows, parse_row(line))
        end
    end
    return rows
end

function time_value(line)
    first_cell = split(line, ';')[1]
    return tryparse(Float64, replace(first_cell, "," => "."))
end

# Read first 3 rows from _1s.csv
rows_1s = read_first_rows(src_path, 3)
row1_1s, row2_1s, row3_1s = rows_1s[1], rows_1s[2], rows_1s[3]
println("_1s.csv columns: $(length(row1_1s))")

# Build lookup: (group, description) -> variable name
lookup = Dict{Tuple{String,String}, String}()
for i in 1:length(row1_1s)
    key = (row1_1s[i], row2_1s[i])
    name = i <= length(row3_1s) ? row3_1s[i] : ""
    if !haskey(lookup, key)
        lookup[key] = name
    end
end
println("Lookup: $(length(lookup)) unique (group, desc) pairs")

# Read all lines from main file
all_lines = readlines(main_path, keep=true)
println("Total lines in main file: $(length(all_lines))")

# Determine whether row 3 is already a names row (non-numeric)
header_lines = all_lines[1:2]
remaining = all_lines[3:end]
has_names_row = isnothing(time_value(remaining[1]))

if has_names_row
    println("Row 3 already has names ('$(split(remaining[1],';')[1])') — keeping it.")
    names_line = [remaining[1]]
    data_start = remaining[2:end]
else
    println("Row 3 is numeric — inserting names row.")
    row1_main = parse_row(all_lines[1])
    row2_main = parse_row(all_lines[2])
    println("main csv columns: $(length(row1_main))")

    function build_new_row(row1, row2, lkp)
        result = String[]
        cnt = 0
        for i in 1:length(row1)
            key = (row1[i], row2[i])
            name = get(lkp, key, "")
            if name != ""; cnt += 1; end
            push!(result, name)
        end
        return result, cnt
    end

    new_row, matched = build_new_row(row1_main, row2_main, lookup)
    println("Matched $matched / $(length(row1_main)) columns with variable names")
    names_line = [join(new_row, ";") * "\n"]
    data_start = remaining
end

# Filter out data rows with t < 0
data_lines = filter(data_start) do line
    t = time_value(line)
    return isnothing(t) || t >= 0.0
end
removed = length(data_start) - length(data_lines)
println("Removed $removed rows with t < 0")

# Write result
new_lines = vcat(header_lines, names_line, data_lines)

open(main_path, "w") do f
    for line in new_lines
        write(f, line)
    end
end
println("Done! File now has $(length(new_lines)) lines.")
