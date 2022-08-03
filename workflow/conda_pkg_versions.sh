# print conda package versions
# only versions of target packages, not downstream dependencies
list=$(mamba list | grep -A99999999 Name)
for pkg in $(grep -m1 -A50000 dependencies "$UKB"/workflow/envs/ukb.env.yml | tail -n +2 | awk '{print $2}') ; do
	echo "$list" | grep "$pkg"
done
