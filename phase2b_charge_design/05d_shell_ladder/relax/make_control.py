import sys
# read a real charges file, output an all-zero version of the SAME length
real=open(sys.argv[1]).read().split(",")
zeros=",".join("0.000" for _ in real)
open(sys.argv[2],"w").write(zeros)
print("wrote %s: %d charges all 0.000 (was %s)"%(sys.argv[2],len(real),sys.argv[1]))
