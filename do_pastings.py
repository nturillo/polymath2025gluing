from paste import *

cases = [14]

for case in cases:
    input_file = f"j4k5/Pj4k5n{case}.txt"
    graphs = read_graphs_from_file(input_file)
    results = get_all_pastes(graphs, min_degree=case-11)
    break

    for size, pastes in results.items():
        output_file = f"pastings/case{case}/d{size}_case{case}_pastes.txt"
        write_graphs_to_file(pastes, output_file)
        print(f"Processed {input_file} for d = {size} and saved pastes to {output_file}.")
