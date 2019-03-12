ls -1 -s -R SGpp | grep -v total > ../new.txt
ls -1 -s -R SGpp_github | grep -v total | sed 's/_github//' > ../old.txt
meld ../old.txt ../new.txt
