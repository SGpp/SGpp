
if __name__ == '__main__':
  filelist =  "../../basefiles.txt"
  f = open(filelist)
  for l in f.readlines():
    code = open(l.strip())
    s = ""
    processed = False
    for codeline in code.readlines():
      if "namespace base" in codeline:
      	processed = True
      if not "namespace sg" in codeline:
        s += codeline
      else:
        if "{" in codeline:
          s += codeline + "namespace base {\n"
        else:
          s += codeline + "{\nnamespace base\n"
      code.close() 
    if not processed:
     codelines = s.split("\n")
     for i in range(len(codelines)-1,-1,-1):
       if "}" in codelines[i]:
         codelines[i] = codelines[i] + "\n}"
         break
     s = "\n".join(codelines)
     codefile = open(l.strip(), "w")
     codefile.write(s)
     codefile.close() 
  f.close()

