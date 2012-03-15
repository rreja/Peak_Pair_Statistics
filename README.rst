If you save your file from Excel, you might have a  character.
To remove  (typed as Control-V, followed by Control M)from your file. use:
$ perl -p -e 's//\n/g;' <your_file> >new_file.txt
