# This is a small script to convert an LMFDB download file into our input data files for totally real fields
# i.e. degN.txt should just be lmfdb_index|coeffs|monogenic, in a compact size (to save on filespace)

# Convert from boolean to smallint (for monogenic)
monogenic_dict = {True : 1, False : -1, None : 0}

f_input = open("lmfdb_nf_fields_0402_2132.sage", 'r')
f_output = open("degN_temp_output.txt", 'w')
f_output.write("lmfdb_index|coeffs|monogenic\n")

for line in f_input:
    if line[0] == '[':
        # Get line of data
        data = line.strip()
        if data[-1]==",":
            data = data[:-1]
        data = eval(data)

        # Get LMFDB index of label
        lmfdb_index = data[0].split('.')[-1]

        # Get polynomial coeffs
        assert(data[1][-1] == 1)
        coeffs = str(data[1][:-1]).replace(' ', '').replace('[', '').replace(']', '')

        # Get whether field is monogenic
        monogenic = str(monogenic_dict[data[2]])

        output_line = lmfdb_index+'|'+coeffs+'|'+monogenic+"\n"
        f_output.write(output_line)

f_input.close()
f_output.close()
