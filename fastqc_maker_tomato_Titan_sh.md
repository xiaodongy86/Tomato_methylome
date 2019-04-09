```{bash}
####
ls synology2/|sed 's#\(.*\)#/usr/local/bin/fastqc -t 4 -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/fastqc /data5/F15FTSUSAT0747_TOMrcwM/synology2/\1/\1_1.fq.gz /data5/F15FTSUSAT0747_TOMrcwM/synology2/\1/\1_2.fq.gz#'
```
