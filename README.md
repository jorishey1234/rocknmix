- Find and write down the right threshold for solid subtraction in the DRY scan: threshold should be between air and solid (th=0.05) .
- Run the python script with the right parameters :
  
>> python3 front.py -S ‘SAMPLE’ -R ‘RESOLUTION’ -T ‘TYPE’ -th 0.05

Other optional parameters are obtain with

>> python3 front.py --help

- Open postproc figures in the sample directory, under ./**** _ post_proc/
- Check if the mean concentration shows a step like (erf function) shape (Front_mean.png)
- Check if concentration histogram show two well separated peaks (Front_hist.png)
- Check if concentration variance shows a single peak that disappear with distance (Front_var.png)

Check figure Front_stats.png

- a Check if the dispersive width grows with downstream distance
- b Check if the mean concentration gradient decrease with distance
- c Check if the max concentration variance decrease with distance
- d Check if the normalized variance is approximately constant
