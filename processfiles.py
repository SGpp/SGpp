import re

def main():
  namespace = "datadriven"
  filelist =  "%sfiles.txt" % namespace
  f = open(filelist)
  for l in f.readlines():
    l = './src/sgpp/' + l.strip()
    code = open(l.strip())
    s = ""
    processed = False
    for codeline in code.readlines():
      if "namespace %s" % namespace in codeline:
        processed = True
      if  re.search("^namespace sg.*", codeline) == None or processed:
        s += codeline
      else:
        if "{" in codeline:
          s += codeline + "namespace %s {\n" % namespace
        else:
          s += codeline + "{\nnamespace %s\n" % namespace
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

if __name__ == '__main__':
	main()


