import glob,os,pickle,sys
sys.set_int_max_str_digits(0)
ch=[]
for p in glob.glob(r"alt2-50000\gd_rigorous_ckpt_50000\c_*.pkl"):
    a,b,_,_=pickle.load(open(p,"rb")); ch.append((int(a),int(b)))
ch.sort()
print("50k chunks:",len(ch),"coverage: [%d, %d)"%(ch[0][0],ch[-1][1]))
print("sizes:",sorted({b-a for a,b in ch})[:5])
print("gaps:",sum(1 for i in range(len(ch)-1) if ch[i][1]!=ch[i+1][0]))
print("first 3:",ch[:3],"last 3:",ch[-3:])
# root b-file count
lines=open("b084945.txt","rb").read().decode("ascii").splitlines()
terms=[l for l in lines if l.strip()]
print("root b-file data lines:",len(terms),"first:",terms[0][:40],"last:",terms[-1][:40])
print("blank tail:",repr(lines[-1]))
