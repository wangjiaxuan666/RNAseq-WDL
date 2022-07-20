# v2.0 2020609

## Change log

1. change the rRNA remove software , abort urmap change to use bowtie2
2. add the new parameter to control the stringtie quantifly whethe use the recreate gtf file.

本次版本更新了核糖体比对软件，之前用的urmap，速度奇慢无比，上次64个15G的样本跑了五天，其中核糖体去除占用了3天
坑人这是！

所以本次更新更换成bowtie2来做核糖体去除.

另外stringtie的定量软件也是坑人的，重构转录本，将基因组已有的基因的定量给安排到自己组装的新基因上，玩我呢。
小鼠多少人的努力构建基因，我做了几十个转录组，一半的基因都被stringtie安排没了。傻逼功能！所以我改了参数regtf，可以选择
是否要用这个重构转录本的方式定量。false就会自动用原始的基因组gtf来定量。

# v2.2 20220627

流程输出的output，会写在软件的位置，其他人一用就报错，修改下，改成在当前运行目录下，创建WDL，并放进去。