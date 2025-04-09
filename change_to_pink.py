input_file = '/Users/pjoglekar/Library/CloudStorage/GoogleDrive-pjoglek@ncsu.edu/My Drive/Huerta_lab/Pseudomonas_Katie_upset/Revised_ANI_dDDH_highlighted.svg'
output_file = '/Users/pjoglekar/Library/CloudStorage/GoogleDrive-pjoglek@ncsu.edu/My Drive/Huerta_lab/Pseudomonas_Katie_upset/Revised_ANI_dDDH_highlighted_pink.svg'

with open(input_file, 'r') as file:
    lines = file.readlines()

with open(output_file, 'w') as file:
    for line in lines:
        if 'bold' in line:
            line = line.replace('#4D4D4D', '#F3359C')
        file.write(line)