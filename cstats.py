import sys,os,re,pickle,json,math
sys.set_int_max_str_digits(0)
def cert(p):
    t=open(p,encoding="ascii").read()
    S=int(re.search(r"S=(\d+)",t).group(1)); lo=int(re.search(r"lam_lo=(\d+)",t).group(1)); hi=int(re.search(r"lam_hi=(\d+)",t).group(1))
    d=[l.strip() for l in t.splitlines() if l.strip().startswith("0.")][-1][2:]
    w=S*math.log10(2)-math.log10(hi-lo)
    return S,lo,hi,len(d),w
for p in [r"40000\golomb_dickman_40000_certified.txt",r"30000\golomb_dickman_30000_certified.txt",r"30000-v2\golomb_dickman_30000_certified.txt",r"30000-v2\golomb_dickman_2000_certified.txt",r"40000\smoke_60.txt"]:
    S,lo,hi,D,w=cert(p)
    print("%-58s S=%5d D=%5d width=10^-%.3f guard=%.3f lo==hi? %s" % (p,S,D,w,w-D,lo==hi))
print("\n30000 vs 30000-v2 endpoints identical:", cert(r"30000\golomb_dickman_30000_certified.txt")[1:3]==cert(r"30000-v2\golomb_dickman_30000_certified.txt")[1:3])
print("\n== meta.pkl 40k ==")
for p in [r"40000\gd_rigorous_ckpt_40000\meta.pkl",r"30000\gd_rigorous_ckpt_30000\meta.pkl"]:
    o=pickle.load(open(p,"rb")); print(p,"->",o)
print("\n== meta_v12 30000-v2 ==")
for p in [r"30000-v2\ckpt_30000\meta_v12.json",r"30000-v2\ckpt_2000\meta_v12.json"]:
    o=json.load(open(p)); print(p,"->",o["parameters"],"chunk",o["chunk"])
print("\n== 50k progress ==")
import glob
n=len(glob.glob(r"alt2-50000\gd_rigorous_ckpt_50000\c_*.pkl"))
KL,KR,ch=149075,158504,100
exp=(KL+KR+1+ch-1)//ch
print("chunks present",n,"expected",exp,"%.1f%%"%(100*n/exp))
