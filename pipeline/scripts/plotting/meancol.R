require(yaml)
require(ggplot2)

avg_cols <- function(c1, c2){
    col = convertColor((c1 + c2) / 2, from="Lab", to="sRGB")
    return(rgb(col))
}

col_scale <- function(){
    x = data.frame(yaml.load_file("colors.yaml")$colors)
    y = col2rgb(as.character(as.matrix(x)))/255

    lab_cols = apply(y, 2, convertColor, from="sRGB", to="Lab", clip=NA)
    z = names(x)
    n = length(z)

    res = c()
    for(i in 1:n){
	for(j in 1:n){
	    s1 = z[i]
	    s2 = z[j]
	    pop = ifelse(s1==s2, s1, sprintf("%s%s", s1, s2))
	    col = avg_cols(lab_cols[,i], lab_cols[,j])
	    v <- c(pop, col)
	    
	    res <- rbind(res, v)
	}
    } 
    rr = res[,2]
    names(rr) = res[,1]

    sc = scale_color_manual(values=rr, aesthetics = c("colour", "fill"))
}
    #saveRDS(sc, "scripts/plotting/cols.rds")
