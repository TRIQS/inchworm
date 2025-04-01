import argparse

def generate_hopping_matrix_file(flavors, t, filename="hopping_matrix.txt"):
    """
    Generate a hopping matrix input file.

    Parameters:
    - flavors: int, total number of flavors (FLAVORS).
    - t: float, hopping parameter.
    - filename: str, output file name.
    """
    with open(filename, "w") as file:
        for i in range(flavors):
            for j in range(flavors):
                # Extract spin labels (i_sigma and j_sigma)
                i_sigma = i % 2  # Spin label for flavor i (0 or 1)
                j_sigma = j % 2  # Spin label for flavor j (0 or 1)
                
                # Set hopping term H_ij
                if i_sigma == j_sigma and i != j:
                    real_part = -1*t
                else:
                    real_part = 0.0
                imag_part = 0.0  # No imaginary part in this case
                
                # Write to file
                file.write(f"{i} {j} {real_part} {imag_part}\n")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a hopping matrix input file.")
    parser.add_argument("--flavors", type=int, required=True, help="Total number of flavors (FLAVORS).")
    parser.add_argument("--t", type=float, required=True, help="Hopping parameter.")
    parser.add_argument("--output", type=str, default="hopping_matrix.txt", help="Output file name.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Generate the hopping matrix file
    generate_hopping_matrix_file(args.flavors, args.t, args.output)
    print(f"Hopping matrix file '{args.output}' generated successfully.")

if __name__ == "__main__":
    main()