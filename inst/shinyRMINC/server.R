
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#

library(shiny)
library(RMINC)
library(plotrix)
library(ggplot2)
library(reshape2)

# theme_black, stolen from: https://jonlefcheck.net/2013/03/11/black-theme-for-ggplot2-2/
# with a few minor modifications to make it look more like theme_classic
theme_black=function(base_size=12,base_family="") {
  theme_grey(base_size=base_size,base_family=base_family) %+replace%
    theme(
      # Specify axis options
      axis.line=element_line(colour="white"), #element_blank(), 
      axis.text.x=element_text(size=base_size*0.8,color="white",
                               lineheight=0.9,vjust=1), 
      axis.text.y=element_text(size=base_size*0.8,color="white",
                               lineheight=0.9,hjust=1), 
      axis.ticks=element_line(color="white",size = 0.2), 
      axis.title.x=element_text(size=base_size,color="white",vjust=1, margin=margin(5,0,0,0)), 
      axis.title.y=element_text(size=base_size,color="white",angle=90,
                                vjust=0.5, margin=margin(0,5,0,0)), 
      axis.ticks.length=unit(0.3,"lines"), 
      axis.ticks.margin=unit(0.5,"lines"),
      # Specify legend options
      legend.background=element_rect(color=NA,fill="black"), 
      legend.key=element_rect(color=NA, fill="black"), 
      legend.key.size=unit(1.2,"lines"), 
      legend.key.height=NULL, 
      legend.key.width=NULL,     
      legend.text=element_text(size=base_size*0.8,color="white"), 
      legend.title=element_text(size=base_size*0.8,face="bold",hjust=0,
                                color="white"), 
      legend.position="right", 
      legend.text.align=NULL, 
      legend.title.align=NULL, 
      legend.direction="vertical", 
      legend.box=NULL,
      # Specify panel options
      panel.background=element_rect(fill="black",color = NA), 
      panel.border=element_blank(), #element_rect(fill=NA,color="white"), 
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank(), 
      panel.margin=unit(0.25,"lines"),  
      # Specify facetting options
      strip.background=element_rect(fill="grey30",color="grey10"), 
      strip.text.x=element_text(size=base_size*0.8,color="white"), 
      strip.text.y=element_text(size=base_size*0.8,color="white",
                                angle=-90), 
      # Specify plot options
      plot.background=element_rect(color="black",fill="black"), 
      plot.title=element_text(size=base_size*1.2,color="white"), 
      plot.margin=unit(c(1,1,0.5,0.5),"lines")
    )
}

#load("preparedData.RData")
# get the data from the calling environment
statsList <- get("statsList", sys.frame(1))
d <- get("d", sys.frame(1))
anatVol <- get("anatVol", sys.frame(1))
gfs <- get("gfs", sys.frame(1))
m <- get("m", sys.frame(1))
modelfunc <- get("modelfunc", sys.frame(1))
vols <- get("volumes", sys.frame(1))
anatLow <- get("anatLow", sys.frame(1))
anatHigh <- get("anatHigh", sys.frame(1))
fdr <- get("fdr", sys.frame(1))
tholds <-
  `if`(is.null(fdr), NULL, thresholds(fdr))
cat(names(statsList))

shinyServer(function(input, output, clientData, session) {

  # update some elements based on the data being accessed
  updateSelectInput(session, "statistic", choices=names(statsList))

  statsChoices <- colnames(gfs)[!( colnames(gfs) %in% c("filenames", "vols"))]
  updateSelectInput(session, "xvar", choices=statsChoices)
  updateSelectInput(session, "colour", choices=statsChoices)
  updateSelectInput(session, "fill", choices=statsChoices)

  # a function to plot either voxels or volumes based on the input data, colours, fill, etc.
  graphData <- function(ydata, ycaption, inputdata) {
    inputdata$data = ydata
    inputdata$xvar <- inputdata[,input$xvar]
    inputdata$colour <- inputdata[,input$colour]
    inputdata$fill <- inputdata[,input$fill]

    #if (statsList[[input$statistic]]$fill == FALSE) {
    #  p <- qplot(xvar, ydata, data=inputdata, geom=input$graphType, colour=colour) +
    #    xlab(statsList[[input$statistic]]$xvar) +
    #    ylab(ycaption) +
    #    scale_colour_brewer(statsList[[input$statistic]]$colour, palette="Set1")
    #}
    #else {
    #  inputdata$fill <- inputdata[,statsList[[input$statistic]]$fill]
    p <- qplot(xvar, ydata, data=inputdata, geom=input$graphType, colour=colour, fill=fill) +
        xlab(input$xvar) +
        ylab(ycaption) +
        scale_colour_brewer(input$colour, palette="Set1") +
        scale_fill_grey(input$fill)
    #}
    return(p)
  }

  getLocation3 <- function() {
    location <- c(input$click_axial$x, input$click_axial$y)
    #location <- ceiling(location * c(d[1], d[2]))
    location <- c(v$loc1, location[2], location[1])
    return(location)
  }

  getLocation2 <- function() {
    location <- c(input$plot_click$x, input$plot_click$y)
    #location <- ceiling(location * c(d[1], d[3]))
    location <- c(location[2], v$loc2,location[1])
    return(location)
  }

  getLocation1 <- function(){
    location <- c(input$click_sagittal$x, input$click_sagittal$y)
    #location <- ceiling(location * c(d[2], d[3]))
    location <- c(location[2],location[1],v$loc3)
    return(location)
  }
  getVoxel <- function() {
    cat("VOXEL", v$loc1, v$loc2, v$loc3, "\n")
    cat("fnames", gfs$filenames[1])
    voxel <- mincGetVoxel(gfs$filenames, v$loc1, v$loc2, v$loc3)
    return(voxel)
  }

  # create a local copy of input data
  rvoxel <- reactive({ getVoxel() })
  rlocation2 <- reactive({ getLocation2() })
  rlocation1 <- reactive({ getLocation1() })
  rlocation3 <- reactive({ getLocation3() })

  v <- reactiveValues(
    loc1 = ceiling(d[3]/2),
    loc2 = ceiling(d[2]/2),
    loc3 = ceiling(d[1]/2),
    voxel = NULL
  )

  observeEvent(input$plot_click, {
    location <- rlocation2()
    v$loc1 <- location[1]
    v$loc2 <- location[2]
    v$loc3 <- location[3]
    if(input$updatePlot) { v$voxel <- rvoxel() }
  })
  observeEvent(input$click_sagittal, {
    location <- getLocation1()
    v$loc1 <- location[1]
    v$loc2 <- location[2]
    v$loc3 <- location[3]
    if (input$updatePlot) { v$voxel <- rvoxel() }
  })
  observeEvent(input$click_axial, {
    location <- rlocation3()
    v$loc1 <- location[1]
    v$loc2 <- location[2]
    v$loc3 <- location[3]
    if (input$updatePlot) { v$voxel <- rvoxel() }
  })

  # update user interface elements based on other selections
  observe({
    # update min and max slice based on the dimension being displaced
    dval <- as.integer(input$dimension)
    diff20 <- round(d[dval] * 0.2)
    cat("in observe slice range:", dval, d[dval], diff20, "\n")
    updateSliderInput(session, "sliceRange", max=d[dval], value=c(diff20, d[dval]-diff20))
    
  })
  observe({
    currentStat <- statsList[[input$statistic]]$data
    currentStat[!is.finite(currentStat)] <- 0 # Maybe not the right thing to do?
    maxstat <- round(max(abs(range(currentStat))), 1)
    low_thresh <-
      as.numeric(
        quantile(abs(currentStat[currentStat != 0])
               , probs = .5))
    
    if (!is.infinite(maxstat)) {
      updateSliderInput(session, "range", max=maxstat, val = c(low_thresh, maxstat))
    }
  })

  output$fdr_thresholds <-
    renderUI({
      if(!is.null(fdr))
        selectInput("fdr_thresh", "Select FDR Limit"
                  , c("none"
                    , paste(rownames(tholds),
                            sprintf("%.3f", tholds[,input$statistic])
                          , sep = " - "))
                   , selected = "none")
    })

    observeEvent(input$fdr_thresh, {
      if(!is.null(fdr) && input$fdr_thresh != "none"){
        thresh <- sub(" - .*", "", input$fdr_thresh)
        updateSliderInput(session, "range", value = c(tholds[thresh, input$statistic], input$range[2]))
      }
    })

  output$seriesPlot <- renderPlot({
    mincPlotSliceSeries(anatVol, statsList[[input$statistic]]$data,
                        anatLow=anatLow, anatHigh=anatHigh, low=input$range[1], 
                        high=input$range[2],
                        begin=input$sliceRange[1], end=input$sliceRange[2], plottitle = input$statistic,
                        dim=as.integer(input$dimension), symmetric=statsList[[input$statistic]]$symmetric,
                        legend=statsList[[input$statistic]]$legendTitle, mfrow=c(input$rows, input$columns))


  }, bg="black")

  output$coronalPlot <- renderPlot({

    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=anatLow, anatHigh=anatHigh,
                              low=input$range[1], high=input$range[2],
                              slice=v$loc2, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=2)
    #mincContour(statsList[[input$statistic]]$data, slice=v$loc2, col="green", add=T, levels=statsList[[input$statistic]]$thresholds[3])
    #abline(h=v$loc1/d[3])
    abline(h=v$loc1)
    #abline(v=v$loc3/d[1])
    abline(v=v$loc3)
  }, bg="black")
  output$sagittalPlot <- renderPlot({
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=anatLow, anatHigh=anatHigh,
                              low=input$range[1], high=input$range[2],
                              slice=v$loc3, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=1, legend="F-statistic", legendTextColour = "white")
    #abline(h=v$loc1/d[3])
    abline(h=v$loc1)
    #abline(v=v$loc2/d[2])
    abline(v=v$loc2)
  }, bg="black")
  output$axialPlot <- renderPlot({
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=anatLow, anatHigh=anatHigh,
                              low=input$range[1], high=input$range[2],
                              slice=v$loc1, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=3)
    #abline(v=v$loc3/d[1])
    abline(v=v$loc3)
    #abline(h=v$loc2/d[2])
    abline(h=v$loc2)
  }, bg="black")

  output$graphPlot <- renderPlot({
    update_geom_defaults("point", list(colour = "white"))
    update_geom_defaults("errorbar", list(colour = "white"))
    update_geom_defaults("boxplot", list(colour = "white", fill="transparent"))
    
    gfs$voxel <- v$voxel
    if(!is.null(gfs$voxel))
      graphData(exp(v$voxel), "jacobians", gfs) + theme_black() + theme(legend.position="top") + 
      guides(colour = guide_legend(nrow=1), fill=guide_legend(nrow=1))

    #qplot(xvar, exp(voxel), data=gfs, geom=input$graphType, colour=colour, fill=fill)
  })

  output$summaryText <- renderPrint({
    #location <- c(input$plot_click$x, input$plot_click$y)
    #location <- ceiling(location * c(d[1], d[3]))
    gfs$voxel <- v$voxel #mincGetVoxel(gfs$filenames, location[2], input$slice, location[1])
    cat("The current coordinates are: ", v$loc1, ",", v$loc2, ",", v$loc3, "\n")
    if(!is.null(gfs$voxel))
      modelfunc(gfs$voxel)
    #anova(lm(voxel ~ mouse.gender + Neonatal, gfs))
    #statsList[[input$statistic]]$modelfunc(gfs)
  })

  output$volumesTable <- DT::renderDataTable({
    usableNames <- input$xvar #colnames(gfs)[!(colnames(gfs) %in% c("filenames", "vols"))]
    #isFalse <- usableNames %in% FALSE
    #usableNames <- usableNames[!isFalse]
    #gfs$newGrouping <- ""
    #for (i in 1:length(usableNames)) {
      #cat("usableNames", i , usableNames[[i]], " ")
      #gfs$newGrouping <- paste(gfs$newGrouping, gfs[,usableNames[[i]]])
    #}
    #cat("done\n")

    #gfs$newGrouping <- paste(gfs[,as.vector(usableNames)], sep="::")
    #gfs$newGrouping <- factor(gfs$newGrouping)
    gfs$newGrouping <- factor(gfs[,usableNames])
    cat(gfs$newGrouping[1])

    dt <- as.data.frame(t(apply(gfs$vols, 2, function(x) { tapply(x, gfs[,"newGrouping"], mean)})))
    dtsd <- as.data.frame(t(apply(gfs$vols, 2, function(x) { tapply(x, gfs[,"newGrouping"], sd)})))
    l <- levels(gfs[,"newGrouping"])
    for (i in 2:length(l)) {
      cat("pre dt\n")
      dt[,paste("effect size:", l[i])] <- (dt[,l[i]] - dt[,l[1]]) / dtsd[,l[1]]
      cat("post dt\n")
    }
    #dt$'effect size' <- (dt[,"male"] - dt[,"female"]) / dtsd[,"male"]

    #colnames(qavs) <- paste("q-value", colnames(qavs))
    #dt <- cbind(qavs, dt)
    dt

  })
  output$volumesPlot <- renderPlot({
    selectedStructures <- input$volumesTable_rows_selected

    usableNames <- colnames(gfs)[!(colnames(gfs) %in% c("filenames", "vols"))]
    #isFalse <- usableNames %in% FALSE
    #usableNames <- usableNames[!isFalse]
    #cat(1)
    v <- as.data.frame(gfs$vols)
    cat("V", nrow(v), ncol(v), "\n")
    for (i in 1:length(usableNames)) {
      v[,usableNames] = gfs[,usableNames]
    }
    #usableNames <- unlist(usableNames)
    #cat(2)
    #v[,statsList[[input$statistic]]$xvar] = gfs[,statsList[[input$statistic]]$xvar]
    #cat(" ", usableNames[1], " ")
    m <- melt(v, id.vars=usableNames) #statsList[[input$statistic]]$xvar)
    cat(3)
    m <- subset(m, variable %in% selectedStructures)
    cat(4)
    for (i in 1:length(names(m))) { cat(names(m)[i], " ")}
    cat("5\n")

    graphData(m$value, "Volume", m) + facet_wrap(~variable, scales="free_y")
    #qplot(m[,statsList[[input$statistic]]$xvar], value, colour=m[,statsList[[input$statistic]]$xvar],
    #      geom=input$graphType, data=m, ylab="Volume (mm3)", xlab=input$statistic) +
    #  scale_colour_brewer(input$statistic, palette="Set1") +
    #  facet_wrap(~variable, scales="free_y")
    })
  #output$volumesPlot2 <- renderPlot({
  #  selectedStructure <- input$volumesTable_rows_selected[length(input$volumesTable_rows_selected)-1]
  #  qplot(mouse.gender, vols[,selectedStructure],
  #        geom="boxplot", data=gfs, main=selectedStructure, ylab="Volume (mm3)")
  #})
})
