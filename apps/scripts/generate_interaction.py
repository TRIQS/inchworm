import argparse

def generate_coulomb_tensor_file(flavors, U, filename="coulomb_tensor.txt"):
    """
    Generate a Coulomb tensor input file.

    Parameters:
    - flavors: int, total number of flavors (FLAVORS).
    - U: float, Coulomb interaction parameter.
    - filename: str, output file name.
    """
    non_zero_elements = []

    # Iterate over all possible combinations of i, j, k, l
    for i in range(flavors):
        for j in range(flavors):
            for k in range(flavors):
                for l in range(flavors):
                    # Extract orbital and spin labels
                    i_orbital = i // 2  # Orbital label for flavor i
                    i_spin = i % 2      # Spin label for flavor i (0 or 1)
                    j_orbital = j // 2  # Orbital label for flavor j
                    j_spin = j % 2      # Spin label for flavor j (0 or 1)
                    k_orbital = k // 2  # Orbital label for flavor k
                    k_spin = k % 2      # Spin label for flavor k (0 or 1)
                    l_orbital = l // 2  # Orbital label for flavor l
                    l_spin = l % 2      # Spin label for flavor l (0 or 1)

                    # Check the conditions:
                    # 1. All indices belong to the same orbital.
                    # 2. i and l have the same spin.
                    # 3. j and k have the opposite spin of i and l.
                    if (i_orbital == j_orbital == k_orbital == l_orbital) and (i_spin == l_spin) and (j_spin == k_spin) and (i_spin != j_spin):
                        non_zero_elements.append((i, j, k, l, U))

    # Write to file
    with open(filename, "w") as file:
        # First line: number of non-zero elements
        file.write(f"{len(non_zero_elements)}\n")
        
        # Subsequent lines: non-zero elements
        for idx, (i, j, k, l, U_val) in enumerate(non_zero_elements):
            file.write(f"{idx} {i} {j} {k} {l} {U_val} 0.0\n")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a Coulomb tensor input file.")
    parser.add_argument("--flavors", type=int, required=True, help="Total number of flavors (FLAVORS).")
    parser.add_argument("--U", type=float, required=True, help="Coulomb interaction parameter.")
    parser.add_argument("--output", type=str, default="coulomb_tensor.txt", help="Output file name.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Generate the Coulomb tensor file
    generate_coulomb_tensor_file(args.flavors, args.U, args.output)
    print(f"Coulomb tensor file '{args.output}' generated successfully.")

if __name__ == "__main__":
    main()