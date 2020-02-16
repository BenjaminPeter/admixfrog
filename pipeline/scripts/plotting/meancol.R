require(yaml)
require(ggplot2)
require(dplyr)

avg_cols <- function(c1, c2){
    col = convertColor((c1 + c2) / 2, from="Lab", to="sRGB")
    return(rgb(col))
}

col_vec <- function(){
    x = data.frame(yaml.load_file("colors.yaml", eval.expr=F)$colors)
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
        if( (pop %in% z) & s1 != s2){next}
	    col = avg_cols(lab_cols[,i], lab_cols[,j])
	    v <- c(pop, col)
	    
	    res <- rbind(res, v)
	}
    } 
    colnames(res) = c("id", "col")
    res = res %>% as_tibble %>% group_by(id) %>% summarize(col=first(col))
    rr = res$col
    names(rr) = res$id
    #print(rr)
    return(rr)

}

col_scale <- function(){
    rr = col_vec()
    sc = scale_color_manual(values=rr, aesthetics = c("colour", "fill"))
}
    #saveRDS(sc, "scripts/plotting/cols.rds")
