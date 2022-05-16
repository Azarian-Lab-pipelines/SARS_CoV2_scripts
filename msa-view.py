def PrintAligns(r,q, text):
    """
    Function to print out the string alignments wtih colors
    """

    r = np.array(list(r))
    q = np.array(list(q))
    
    if not text:
        print(f"Ref: {''.join(r)}")
        print("Qry: " + "".join([f"{bcolors.MIS}{q[n]}{bcolors.RESET}" if _ == False \
            else q[n]
            for n,_ in enumerate(q == r)]))
        print()
    else:
        print("Ref: " + "".join([f"{r[n]}" if _ == False \
            else "*"
            for n,_ in enumerate(q == r)]))
        print("Qry: " + "".join([f"{q[n]}" if _ == False \
            else "*"
            for n,_ in enumerate(q == r)]))
        print()

# Just figure out how many characters per line to print
F = 
