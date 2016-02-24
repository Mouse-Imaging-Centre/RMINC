
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

#load("preparedData.RData")
# get the data from the calling environment
statsList <- get("statsList", sys.frame(1))
d <- get("d", sys.frame(1))
anatVol <- get("anatVol", sys.frame(1))
gfs <- get("gfs", sys.frame(1))
m <- get("m", sys.frame(1))
modelfunc <- get("modelfunc", sys.frame(1))
vols <- get("volumes", sys.frame(1))
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
    updateSliderInput(session, "begin", max=d[dval])
    updateSliderInput(session, "end", min=-d[dval])
  })
  observe({
    currentStat <- input$statistic
    maxstat <- round(max(abs(range(statsList[[currentStat]]$data))), 1)
    cat("in observe ", maxstat, "\n")
    if (!is.infinite(maxstat)) {
      updateSliderInput(session, "high", max=maxstat)
      updateSliderInput(session, "low", max=maxstat)
    }
  })

  output$seriesPlot <- renderPlot({
    cat("Low", input$low, "High", input$high, "sym", statsList[[input$statistic]]$symmetric, "\n")
    mincPlotSliceSeries(anatVol, statsList[[input$statistic]]$data,
                        anatLow=700, anatHigh=1400, low=input$low, high=input$high,
                        begin=input$begin, end=input$end, plottitle = input$statistic,
                        dim=as.integer(input$dimension), symmetric=statsList[[input$statistic]]$symmetric,
                        legend=statsList[[input$statistic]]$legendTitle, mfrow=c(input$rows, input$columns))


  }, bg="black")

  output$coronalPlot <- renderPlot({

    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=700, anatHigh=1400,
                              low=input$low, high=input$high,
                              slice=v$loc2, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=2)
    #mincContour(statsList[[input$statistic]]$data, slice=v$loc2, col="green", add=T, levels=statsList[[input$statistic]]$thresholds[3])
    #abline(h=v$loc1/d[3])
    abline(h=v$loc1)
    #abline(v=v$loc3/d[1])
    abline(v=v$loc3)
  })
  output$sagittalPlot <- renderPlot({
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=700, anatHigh=1400,
                              low=input$low, high=input$high,
                              slice=v$loc3, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=1, legend="F-statistic")
    #abline(h=v$loc1/d[3])
    abline(h=v$loc1)
    #abline(v=v$loc2/d[2])
    abline(v=v$loc2)
  })
  output$axialPlot <- renderPlot({
    mincPlotAnatAndStatsSlice(anatVol,
                              statsList[[input$statistic]]$data,
                              anatLow=700, anatHigh=1400,
                              low=input$low, high=input$high,
                              slice=v$loc1, symmetric=statsList[[input$statistic]]$symmetric,
                              dim=3)
    #abline(v=v$loc3/d[1])
    abline(v=v$loc3)
    #abline(h=v$loc2/d[2])
    abline(h=v$loc2)
  })

  output$graphPlot <- renderPlot({
    gfs$voxel <- v$voxel
    graphData(exp(v$voxel), "jacobians", gfs)

    #qplot(xvar, exp(voxel), data=gfs, geom=input$graphType, colour=colour, fill=fill)
  })

  output$summaryText <- renderPrint({
    #location <- c(input$plot_click$x, input$plot_click$y)
    #location <- ceiling(location * c(d[1], d[3]))
    gfs$voxel <- v$voxel #mincGetVoxel(gfs$filenames, location[2], input$slice, location[1])
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
