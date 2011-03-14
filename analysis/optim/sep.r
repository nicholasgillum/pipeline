library(stats)

dataframe <- read.table("exp_data.txt")
rawtime   <- dataframe$time
rawdata   <- dataframe$VNG1182H

premask <- rawtime <= 0;

ptime <- rawtime[premask];
pexpr <- rawdata[premask];

time  <- rawtime[!premask];
expr  <- rawdata[!premask];



plot(time, expr)
points(ptime, pexpr)


