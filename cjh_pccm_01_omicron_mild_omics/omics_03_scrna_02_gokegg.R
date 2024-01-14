library(cowplot)
library(grid)
library(gridExtra)

MYCOLOR <- c("#FF6700", "#9370DB", "#F8D568",
            "#00AD43", "#89CFF0", "#BA160C",
            "#FF91AF", "#A6A6A6", "#006DB0",
            "#C154C1", "#D99A6C", "#96C8A2",
            "#FBEC5D", "#77B5FE", "#E29CD2",
            "#007F5C", "#ACBF60", "#7B68EE",
            "#00FFCD", "#D3AF37", "#50C878",
            "#1974D2", "#FF0000", "#9F00FF",
            "#91A3B0", "#8B4513", "#4166F5",
            "#C19A6B", "#6EAEA1", "#39FF14",
            "#0247FE", "#AF6E4D", "#FF66CC",
            "#9400D3", "#A57164", "#00FFFF",
            "#DF00FF", "#FFAA1D", "#20B2AA",
            "#CB6D51", "#5218FA", "#BC8F8F",
            "#000000")
MYCOLOR <- MYCOLOR[-14]



refineDoHeatmap <- function (object, features = NULL, cells = NULL, group.by = "ident",
    group.bar = TRUE, disp.min = -2.5, disp.max = NULL, slot = "scale.data",
    assay = NULL, label = TRUE, size = 5.5, hjust = 0, angle = 45,
    raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02,
    combine = TRUE)
{
    cells <- cells %||% colnames(x = object)
    if (is.numeric(x = cells)) {
        cells <- colnames(x = object)[cells]
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    features <- features %||% VariableFeatures(object = object)
    features <- rev(x = unique(x = features))
    disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
        yes = 2.5, no = 6)
    possible.features <- rownames(x = GetAssayData(object = object,
        slot = slot))
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        features <- features[features %in% possible.features]
        if (length(x = features) == 0) {
            stop("No requested features found in the ", slot,
                " slot for the ", assay, " assay.")
        }
        warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                collapse = ", "))
    }
    data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object,
        slot = slot)[features, cells, drop = FALSE])))
    object <- suppressMessages(expr = StashIdent(object = object,
        save.name = "ident"))
    group.by <- group.by %||% "ident"
    groups.use <- object[[group.by]][cells, , drop = FALSE]
    plots <- vector(mode = "list", length = ncol(x = groups.use))
    for (i in 1:ncol(x = groups.use)) {
        data.group <- data
        group.use <- groups.use[, i, drop = TRUE]
        if (!is.factor(x = group.use)) {
            group.use <- factor(x = group.use)
        }
        names(x = group.use) <- cells
        if (draw.lines) {
            lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                0.0025)
            placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) *
                lines.width), FUN = function(x) {
                return(Seurat:::RandomName(length = 20))
            })
            placeholder.groups <- rep(x = levels(x = group.use),
                times = lines.width)
            group.levels <- levels(x = group.use)
            names(x = placeholder.groups) <- placeholder.cells
            group.use <- as.vector(x = group.use)
            names(x = group.use) <- cells
            group.use <- factor(x = c(group.use, placeholder.groups),
                levels = group.levels)
            na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                  colnames(x = data.group)))
            data.group <- rbind(data.group, na.data.group)
        }
        # internal function SingleRasterMap
        plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster,
            disp.min = disp.min, disp.max = disp.max, feature.order = features,
            cell.order = names(x = sort(x = group.use)), group.by = group.use)

        if (group.bar) {
            group.use2 <- sort(x = group.use)
            if (draw.lines) {
                na.group <- Seurat:::RandomName(length = 20)
                levels(x = group.use2) <- c(levels(x = group.use2),
                  na.group)
                group.use2[placeholder.cells] <- na.group
                cols <- c(MYCOLOR[1:length(x = levels(x = group.use))],
                  "#FFFFFF")
            }
            else {
                cols <- c(MYCOLOR[1:length(x = levels(x = group.use))])
            }
            plot <- plot + geom_point(mapping = aes_string(x = "Cell",
                        y = "Feature", color = "Identity"), alpha = 0) +
                        scale_colour_manual(values = cols)

            pbuild <- ggplot_build(plot = plot)
            names(x = cols) <- levels(x = group.use2)
            y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
            y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) +
                y.range * 0.015
            y.max <- y.pos + group.bar.height * y.range
            plot <- plot + annotation_raster(raster = t(x = cols[group.use2]),
                , xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
                coord_cartesian(ylim = c(0, y.max), clip = "off")
            if (label) {
                x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
                x.divs <- pbuild$layout$panel_params[[1]]$x.major
                x <- data.frame(group = sort(x = group.use),
                  x = x.divs)
                label.x.pos <- tapply(X = x$x, INDEX = x$group,
                  FUN = median) * x.max
                label.x.pos <- data.frame(group = names(x = label.x.pos),
                  label.x.pos)
                plot <- plot + geom_text(stat = "identity", data = label.x.pos,
                  aes_string(label = "group", x = "label.x.pos"),
                  y = y.max + y.max * 0.03 * 0.5, angle = angle,
                  hjust = hjust, size = size)
                plot <- suppressMessages(plot + coord_cartesian(ylim = c(0,
                  y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) *
                    size), clip = "off"))
            }
        }
        plot <- plot + theme(line = element_blank(), axis.text.y =  element_text(size=5))
        plots[[i]] <- plot
    }
    if (combine) {
        plots <- CombinePlots(plots = plots)
    }
    return(plots)
}


arrangeTop9 <- function(plotlist,  legend.position='right') {
    plotlist <- plotlist %>% map(~ { .x +  theme(axis.title=element_blank(), axis.text=element_text(size=rel(0.5)), legend.position=legend.position)})
    p0 <- plot_grid(plotlist=plotlist)
    
    pbuild <- ggplot_build(plot = plotlist[[1]])
    plot_labels_y <- pbuild$plot$labels$y
    plot_labels_x <- pbuild$plot$labels$x

    y.grob <- textGrob(plot_labels_y,
                       gp=gpar(family='sans', fontsize=15), rot=90)
    
    x.grob <- textGrob(plot_labels_x,
                       gp=gpar(family='sans', fontsize=15))
    
    p1 <-grid.arrange(arrangeGrob(p0, bottom=x.grob, left=y.grob))
    return(p1)
}

plotBar <- function(df, n=10, pname='p.adjust'){
    df <- df %>% dplyr::rename("Ontology"="ONTOLOGY") %>% 
              group_by(Ontology) %>% 
              top_n(-n, wt=!!sym(pname)) %>%
              mutate(logp=-(log10(!!sym(pname)))) %>%
              arrange(Ontology, desc(logp))

    df$Description <- factor(df$Description, levels = df$Description)
    p <- ggplot(df, aes(Description, logp, fill=Ontology))
    p <- p + geom_bar(stat="identity") + labs(y=sprintf("-log10(%s)", pname)) + theme_classic() + 
                 theme(plot.margin = margin(5.5, 5.5, 5.5, 100, "pt"), axis.text = element_text(angle=45, hjust=1))
    return(p)
}

plotBar2 <- function(df, n=20, pname='p.adjust'){
    df <- df %>% top_n(-n, wt=!!sym(pname)) %>%
              mutate(logp=-(log10(!!sym(pname)))) %>%
              rowwise() %>%
              mutate(Description2=str_wrap(Description, 40))

    df$Description2 <- factor(df$Description2, levels = df$Description2)
    p <- ggplot(df, aes(Description2, Count, fill=logp))
    p <- p + geom_bar(stat="identity") + labs(x="Description") +  theme_classic() + 
              scale_fill_gradient(sprintf("-log10(%s)", pname), low="blue", high="red") + coord_flip()
    return(p)
}


plotDot <- function(df, n=20, pname='p.adjust'){
    df <- df %>% top_n(-n, wt=!!sym(pname)) %>% rowwise() %>% mutate(GeneRatio2=eval(parse(text=GeneRatio))) %>% 
              mutate(Description2=str_wrap(Description, 40)) %>%
              arrange(desc(GeneRatio2))
    df$Description2 <- factor(df$Description2, levels=df$Description2)
    p <- ggplot(df, aes(GeneRatio2, Description2, color=!!sym(pname), size=Count))
    p <- p + geom_point() + scale_colour_gradient(pname, limits=c(0, 1), high="blue", low="red") +
                 labs(x="GeneRatio", y="Description") +
                 guides(size = guide_legend(order=1), colour = guide_colorbar(order=2))
    return(p)
}

pathview.new <- function (gene.data = NULL, cpd.data = NULL, pathway.id, species = "hsa",
    kegg.dir = ".", cpd.idtype = "kegg", gene.idtype = "entrez",
    gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE,
    map.null = TRUE, expand.node = FALSE, split.group = FALSE,
    map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum",
    discrete = list(gene = FALSE, cpd = FALSE), limit = list(gene = 1,
        cpd = 1), bins = list(gene = 10, cpd = 10), both.dirs = list(gene = T,
        cpd = T), trans.fun = list(gene = NULL, cpd = NULL),
    low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray",
        cpd = "gray"), high = list(gene = "red", cpd = "yellow"),
    na.col = "transparent", ...)
{
    dtypes = !is.null(gene.data) + !is.null(cpd.data)
    cond0 = dtypes == 1 & is.numeric(limit) & length(limit) >
        1
    if (cond0) {
        if (limit[1] != limit[2] & is.null(names(limit)))
            limit = list(gene = limit[1:2], cpd = limit[1:2])
    }
    if (is.null(trans.fun))
        trans.fun = list(gene = NULL, cpd = NULL)
    arg.len2 = c("discrete", "limit", "bins", "both.dirs", "trans.fun",
        "low", "mid", "high")
    for (arg in arg.len2) {
        obj1 = eval(as.name(arg))
        if (length(obj1) == 1)
            obj1 = rep(obj1, 2)
        if (length(obj1) > 2)
            obj1 = obj1[1:2]
        obj1 = as.list(obj1)
        ns = names(obj1)
        if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns))
            names(obj1) = c("gene", "cpd")
        assign(arg, obj1)
    }
    if (is.character(gene.data)) {
        gd.names = gene.data
        gene.data = rep(1, length(gene.data))
        names(gene.data) = gd.names
        both.dirs$gene = FALSE
        ng = length(gene.data)
        nsamp.g = 1
    }
    else if (!is.null(gene.data)) {
        if (length(dim(gene.data)) == 2) {
            gd.names = rownames(gene.data)
            ng = nrow(gene.data)
            nsamp.g = 2
        }
        else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
            gd.names = names(gene.data)
            ng = length(gene.data)
            nsamp.g = 1
        }
        else stop("wrong gene.data format!")
    }
    else if (is.null(cpd.data)) {
        stop("gene.data and cpd.data are both NULL!")
    }
    gene.idtype = toupper(gene.idtype)
    data(bods)
    if (species != "ko") {
        species.data = kegg.species.code(species, na.rm = T,
            code.only = FALSE)
    }
    else {
        species.data = c(kegg.code = "ko", entrez.gnodes = "0",
            kegg.geneid = "K01488", ncbi.geneid = NA, ncbi.proteinid = NA,
            uniprot = NA)
        gene.idtype = "KEGG"
        msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
        msg = sprintf(msg.fmt, species.data["kegg.geneid"])
        message("Note: ", msg)
    }
    if (length(dim(species.data)) == 2) {
        message("Note: ", "More than two valide species!")
        species.data = species.data[1, ]
    }
    species = species.data["kegg.code"]
    entrez.gnodes = species.data["entrez.gnodes"] == 1
    if (is.na(species.data["ncbi.geneid"])) {
        if (!is.na(species.data["kegg.geneid"])) {
            msg.fmt = "Mapping via KEGG gene ID (not Entrez) is supported for this species,\nit looks like \"%s\"!"
            msg = sprintf(msg.fmt, species.data["kegg.geneid"])
            message("Note: ", msg)
        }
        else {
            stop("This species is not annotated in KEGG!")
        }
    }
    if (is.null(gene.annotpkg))
        gene.annotpkg = bods[match(species, bods[, 3]), 1]
    if (length(grep("ENTREZ|KEGG|NCBIPROT|UNIPROT", gene.idtype)) <
        1 & !is.null(gene.data)) {
        if (is.na(gene.annotpkg))
            stop("No proper gene annotation package available!")
        if (!gene.idtype %in% gene.idtype.bods[[species]])
            stop("Wrong input gene ID type!")
        gene.idmap = id2eg(gd.names, category = gene.idtype,
            pkg.name = gene.annotpkg, unique.map = F)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "ENTREZ"
    }
    if (gene.idtype != "KEGG" & !entrez.gnodes & !is.null(gene.data)) {
        id.type = gene.idtype
        if (id.type == "ENTREZ")
            id.type = "ENTREZID"
        kid.map = names(species.data)[-c(1:2)]
        kid.types = names(kid.map) = c("KEGG", "ENTREZID", "NCBIPROT",
            "UNIPROT")
        kid.map2 = gsub("[.]", "-", kid.map)
        kid.map2["UNIPROT"] = "up"
        if (is.na(kid.map[id.type]))
            stop("Wrong input gene ID type for the species!")
        message("Info: Getting gene ID data from KEGG...")
        gene.idmap = keggConv(kid.map2[id.type], species)
        message("Info: Done with data retrieval!")
        kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
        in.ids = gsub(paste0(kid.map2[id.type], ":"), "", gene.idmap)
        gene.idmap = cbind(in.ids, kegg.ids)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "KEGG"
    }
    if (is.character(cpd.data)) {
        cpdd.names = cpd.data
        cpd.data = rep(1, length(cpd.data))
        names(cpd.data) = cpdd.names
        both.dirs$cpd = FALSE
        ncpd = length(cpd.data)
    }
    else if (!is.null(cpd.data)) {
        if (length(dim(cpd.data)) == 2) {
            cpdd.names = rownames(cpd.data)
            ncpd = nrow(cpd.data)
        }
        else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
            cpdd.names = names(cpd.data)
            ncpd = length(cpd.data)
        }
        else stop("wrong cpd.data format!")
    }
    if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
        data(rn.list)
        cpd.types = c(names(rn.list), "name")
        cpd.types = tolower(cpd.types)
        cpd.types = cpd.types[-grep("kegg", cpd.types)]
        if (!tolower(cpd.idtype) %in% cpd.types)
            stop("Wrong input cpd ID type!")
        cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
        cpd.data = mol.sum(cpd.data, cpd.idmap)
    }
    warn.fmt = "Parsing %s file failed, please check the file!"
    if (length(grep(species, pathway.id)) > 0) {
        pathway.name = pathway.id
        pathway.id = gsub(species, "", pathway.id)
    }
    else pathway.name = paste(species, pathway.id, sep = "")
    kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
    npath = length(pathway.id)
    out.list = list()
    tfiles.xml = paste(pathway.name, "xml", sep = ".")
    tfiles.png = paste(pathway.name, "png", sep = ".")
    if (kegg.native)
        ttype = c("xml", "png")
    else ttype = "xml"
    xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
    for (i in 1:npath) {
        if (kegg.native)
            tfiles = c(tfiles.xml[i], tfiles.png[i])
        else tfiles = tfiles.xml[i]
        if (!all(tfiles %in% kfiles)) {
            dstatus = download.kegg(pathway.id = pathway.id[i],
                species = species, kegg.dir = kegg.dir, file.type = ttype)
            if (dstatus == "failed") {
                warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
                warn.msg = sprintf(warn.fmt, pathway.name[i])
                message("Warning: ", warn.msg)
                return(invisible(0))
            }
        }
        if (kegg.native) {
            node.data = try(node.info(xml.file[i]), silent = T)
            if (class(node.data)[1] == "try-error") {
                warn.msg = sprintf(warn.fmt, xml.file[i])
                message("Warning: ", warn.msg)
                return(invisible(0))
            }
            node.type = c("gene", "enzyme", "compound", "ortholog")
            sel.idx = node.data$type %in% node.type
            nna.idx = !is.na(node.data$x + node.data$y + node.data$width +
                node.data$height)
            sel.idx = sel.idx & nna.idx
            if (sum(sel.idx) < min.nnodes) {
                warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
                warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
                message("Warning: ", warn.msg)
                return(invisible(0))
            }
            node.data = lapply(node.data, "[", sel.idx)
        }
        else {
            gR1 = try(parseKGML2Graph2(xml.file[i], genes = F,
                expand = expand.node, split.group = split.group),
                silent = T)
            node.data = try(node.info(gR1), silent = T)
            if (class(node.data)[1] == "try-error") {
                warn.msg = sprintf(warn.fmt, xml.file[i])
                message("Warning: ", warn.msg)
                return(invisible(0))
            }
        }
        if (species == "ko")
            gene.node.type = "ortholog"
        else gene.node.type = "gene"
        if ((!is.null(gene.data) | map.null) & sum(node.data$type ==
            gene.node.type) > 1) {
            plot.data.gene = node.map(gene.data, node.data, node.types = gene.node.type,
                node.sum = node.sum, entrez.gnodes = entrez.gnodes)
			plot.data.gene<-plot.data.gene[rowSums(plot.data.gene[,c("x","y","width","height")])!=4,]
            kng = plot.data.gene$kegg.names
            kng.char = gsub("[0-9]", "", unlist(kng))
            if (any(kng.char > ""))
                entrez.gnodes = FALSE
            if (map.symbol & species != "ko" & entrez.gnodes) {
                if (is.na(gene.annotpkg)) {
                  warn.fmt = "No annotation package for the species %s, gene symbols not mapped!"
                  warn.msg = sprintf(warn.fmt, species)
                  message("Warning: ", warn.msg)
                }
                else {
                  plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names),
                    category = "SYMBOL", pkg.name = gene.annotpkg)[,
                    2]
                  mapped.gnodes = rownames(plot.data.gene)
                  node.data$labels[mapped.gnodes] = plot.data.gene$labels
                }
            }
            cols.ts.gene = node.color(plot.data.gene, limit$gene,
                bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene,
                discrete = discrete$gene, low = low$gene, mid = mid$gene,
                high = high$gene, na.col = na.col)
        }
        else plot.data.gene = cols.ts.gene = NULL
        if ((!is.null(cpd.data) | map.null) & sum(node.data$type ==
            "compound") > 1) {
            plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound",
                node.sum = node.sum)
            if (map.cpdname & !kegg.native) {
                plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[,
                  2]
                mapped.cnodes = rownames(plot.data.cpd)
                node.data$labels[mapped.cnodes] = plot.data.cpd$labels
            }
            cols.ts.cpd = node.color(plot.data.cpd, limit$cpd,
                bins$cpd, both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd,
                discrete = discrete$cpd, low = low$cpd, mid = mid$cpd,
                high = high$cpd, na.col = na.col)
        }
        else plot.data.cpd = cols.ts.cpd = NULL
        if (kegg.native) {
            pv.pars = keggview.native(plot.data.gene = plot.data.gene,
                cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd,
                cols.ts.cpd = cols.ts.cpd, node.data = node.data,
                pathway.name = pathway.name[i], kegg.dir = kegg.dir,
                limit = limit, bins = bins, both.dirs = both.dirs,
                discrete = discrete, low = low, mid = mid, high = high,
                na.col = na.col, ...)
        }
        else {
            pv.pars = keggview.graph(plot.data.gene = plot.data.gene,
                cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd,
                cols.ts.cpd = cols.ts.cpd, node.data = node.data,
                path.graph = gR1, pathway.name = pathway.name[i],
                map.cpdname = map.cpdname, split.group = split.group,
                limit = limit, bins = bins, both.dirs = both.dirs,
                discrete = discrete, low = low, mid = mid, high = high,
                na.col = na.col, ...)
        }
        plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
        if (!is.null(plot.data.gene)) {
            cnames = colnames(plot.data.gene)[-(1:8)]
            nsamp = length(cnames)/2
            if (nsamp > 1) {
                cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp +
                  1):(2 * nsamp)], "col", sep = ".")
            }
            else cnames[2] = "mol.col"
            colnames(plot.data.gene)[-(1:8)] = cnames
        }
        plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
        if (!is.null(plot.data.cpd)) {
            cnames = colnames(plot.data.cpd)[-(1:8)]
            nsamp = length(cnames)/2
            if (nsamp > 1) {
                cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp +
                  1):(2 * nsamp)], "col", sep = ".")
            }
            else cnames[2] = "mol.col"
            colnames(plot.data.cpd)[-(1:8)] = cnames
        }
        out.list[[i]] = list(plot.data.gene = plot.data.gene,
            plot.data.cpd = plot.data.cpd)
    }
    if (npath == 1)
        out.list = out.list[[1]]
    else names(out.list) = pathway.name
    return(invisible(out.list))
}
