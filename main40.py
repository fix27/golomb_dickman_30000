import re
p=r"40000\golomb_dickman.py"
t=open(p,encoding="utf-8",errors="replace").read()
i=t.find("def main()")
print("=== 40000 main ===")
print(t[i:i+1200])
print("\n=== run_certificate signature/certify calls ===")
for m in re.finditer(r"run_certificate\(",t):
    print(t[m.start():m.start()+220].replace("\n"," ")[:220])
