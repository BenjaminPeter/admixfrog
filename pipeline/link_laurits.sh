FOLDER=/r1/people/laurits_skov/Admixfrog/infiles/

for file in `ls $FOLDER`; do
    fname=`echo $file | cut -f1 -d.`
    newname="${fname}_archaicadmixture.in.xz"
    echo $fname $newname
    echo $fname >> config/panels.yaml
    ln -s $FOLDER/$file samples/$newname
done

