from string import Template

def main():
  f = open('run_job.tpl', 'r')
  try:

    tpl = Template(f.read())
    for cond in [0,1]:
      for num in [1,5,10]:
        for ref_type in [0,1,2,3]:
          job = tpl.safe_substitute(cond=cond, num=num, ref_type=ref_type)
          o = open("jobs/cond%d-num%d-ref_type%d.job"%(cond, num, ref_type), 'wt')
          o.write(job)
          o.close()
#          break
#        break
#      break

  finally:
    f.close()

if __name__ == "__main__":
	main()

