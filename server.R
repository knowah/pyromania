library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(reshape2)
source("pyro_data_tools.R")

function(input, output, session) {
  ready.for.analysis <- reactiveVal(FALSE)
  PLATE.ROWS <- 8
  PLATE.COLS <- 12

  plate.df <- data.frame(r=rep(1:PLATE.ROWS, each=PLATE.COLS), c=1:PLATE.COLS)
  plate.df$w <- with(plate.df, paste0(LETTERS[r], c))
  plate.df$s <- ""
  
  r.plate <- reactiveVal(plate.df)
  r.pyro.long <- reactiveVal()
  r.assay.summary <- reactiveVal()
  r.CG.select <- reactiveVal()
  r.input.file <- reactiveVal()
  
  # number of samples changed
  observeEvent(input$N_samps, {
    n.samps <- input$N_samps
    if (!is.null(n.samps) && !is.na(n.samps)) {
      n.samps <- as.integer(n.samps)
      if (n.samps > 0) {
        updateSelectInput(session, "samp_select",
                          choices = c(1:n.samps, "Water")
        )
      }
    }
  })
  
  # 'Start Analysis' button clicked
  observeEvent(input$start_analysis, {
    if ("Assay" %in% colnames(r.plate())) {
      ready.for.analysis(TRUE)
      pyro.long <- melt(
        r.plate() %>%
          filter(nchar(as.character(s)) > 0 & nchar(as.character(Assay)) > 0) %>%
          dplyr::select(w, s, Assay, ends_with("Methylation")),
        id.vars = c("w", "s", "Assay"),
        variable.name = "pos",
        value.name = "methy"
      )
      
      pyro.long$pos <- gsub("Pos\\.([0-9]+)\\.Methylation", "\\1", pyro.long$pos)
      
      pyro.long.q <- melt(
        r.plate() %>%
          filter(nchar(as.character(s)) > 0 & nchar(as.character(Assay)) > 0) %>%
          dplyr::select(w, s, Assay, ends_with("Quality")),
        id.vars = c("w", "s", "Assay"),
        variable.name = "pos",
        value.name = "quality"
      )
      
      pyro.long.q$pos <- gsub("Pos\\.([0-9]+)\\.Quality", "\\1", pyro.long.q$pos)
      pyro.long <- left_join(pyro.long, pyro.long.q, by=c("w", "s", "Assay", "pos"))
      
      # remove entries corresponding to unanalyzed CpGs (i.e. CpGs not present on an assay)
      unanalyzed.idx <- with(pyro.long, which(is.na(methy) & nchar(quality) == 0))
      if (length(unanalyzed.idx) > 0) {
        pyro.long <- pyro.long[-unanalyzed.idx,]
      }
      
      pyro.long$exclude <- FALSE
      
      r.pyro.long(pyro.long)
      
      active.assays <- unique(as.character((r.plate() %>% filter(nchar(as.character(s)) > 0) %>% arrange(r, c))$Assay))
      
      max.CG <- 1
      valid.CG <- if (is.null(r.CG.select())) list() else r.CG.select()
      for (assay in active.assays) {
        cgs <- unique((pyro.long %>% filter(Assay==assay & !is.na(methy)))$pos)
        if (!assay %in% names(valid.CG)) {
          valid.CG[[assay]] <- min(cgs):max(cgs)
        }
        max.CG <- max(max.CG, max(cgs))
      }
      
      updateCheckboxGroupInput(session, "CpG_select", choices = 1:max.CG)
      r.CG.select(valid.CG)
      updateSelectInput(session, "assay_select", choices=active.assays, selected = active.assays[1])
      updateCheckboxGroupInput(session, "CpG_select", selected = r.CG.select()[[input$assay_select]])
      
    } 
  })
  
  get.pyro.summary <- function() {
    req(r.pyro.long())
    
    # summarize pyro data, excluding marked samples
    assay.summary <- r.pyro.long() %>%
      filter(!exclude) %>%
      group_by(s, Assay, pos) %>%
      summarize(m.mean=mean(methy, na.rm=TRUE), m.sd=sd(methy, na.rm=TRUE), m.N=sum(!is.na(methy))) %>%
      filter(m.N > 0)
    
    # exclude CpGs not wanted
    excluded.CG.idx <- unlist(sapply(names(r.CG.select()),
                                     function(n) with(assay.summary, which(Assay == n & ! pos %in% r.CG.select()[[n]]))))
    if (length(excluded.CG.idx) > 0)
      assay.summary <- assay.summary[-excluded.CG.idx,]
    
    assay.summary
  }
  
  # update the current summarized methylation data
  update_summary <- function() {
    #req(r.pyro.long())
    #assay.summary <- r.pyro.long() %>%
    #  filter(!exclude)
    #if (input$remove_failed) {
    #  assay.summary <- assay.summary %>% filter(quality != "Failed")
    #}
    #assay.summary <- assay.summary %>%
    #  filter(Assay == input$assay_select) %>%
    #  filter(pos %in% r.CG.select()[[input$assay_select]]) %>%
    #  group_by(s, pos) %>%
    #  summarize(m.mean=mean(methy, na.rm=TRUE), m.sd=sd(methy, na.rm=TRUE), m.N=sum(!is.na(methy))) %>%
    #  filter(m.N > 0)
    r.assay.summary(get.pyro.summary() %>% filter(Assay == input$assay_select) %>% ungroup() %>% dplyr::select(-Assay))
  }
  
  observeEvent(input$remove_failed, {
    req(r.pyro.long())
    if (any(r.pyro.long()$quality == "Failed")) {
      pyro.long <- r.pyro.long()
      pyro.long[pyro.long$quality == "Failed",]$exclude <- TRUE
      r.pyro.long(pyro.long)
      update_summary()
    }
  })
  
  # CpGs to display altered
  observeEvent(input$CpG_select, {
    cg_s <- r.CG.select()
    cg_s[[input$assay_select]] <- input$CpG_select
    r.CG.select(cg_s)
    update_summary()
  })
  
  # Assay to display altered
  observeEvent(input$assay_select, {
    req(r.pyro.long())
    updateCheckboxGroupInput(session, "CpG_select", selected=r.CG.select()[[input$assay_select]])
    update_summary()
  })
  
  # New data file uploaded
  observeEvent(input$pyro_file, {
    plate.df <- r.plate()
    r.input.file(input$pyro_file$name)
    pyro.data <- read.pyro.data(input$pyro_file$datapath)
    r.plate(
      left_join(plate.df, pyro.data, by=c("w"="Well"))
    )
  })
  
  assign.sample <- function(df.rows.to.change, selected.sample) {
    plate.df <- r.plate()
    if (length(df.rows.to.change) > 0) {
      existing.samples <- unique(plate.df[df.rows.to.change,]$s) 
      if (length(existing.samples) == 1 && existing.samples == selected.sample) {
        plate.df[df.rows.to.change,]$s <- ""
      } else {
        plate.df[df.rows.to.change,]$s <- selected.sample 
      }
      r.plate(plate.df)
    }
  }
  
  # Plate wells selected
  observeEvent(input$plate_brush, {
    plate.df <- r.plate()
    
    current.sample <- input$samp_select
    
    r0 <- max(1, floor(input$plate_brush$ymin))
    r1 <- min(PLATE.ROWS, floor(input$plate_brush$ymax))
    c0 <- max(1, floor(input$plate_brush$xmin))
    c1 <- min(PLATE.COLS, floor(input$plate_brush$xmax))
    
    df.rows.to.change <- with(plate.df, which(r >= r0 & r <= r1 & c >= c0 & c <= c1))
    
    if (length(df.rows.to.change) > 0) {
      existing.samples <- unique(plate.df[df.rows.to.change,]$s) 
      if (length(existing.samples) == 1 && existing.samples == current.sample) {
        plate.df[df.rows.to.change,]$s <- ""
      } else {
        plate.df[df.rows.to.change,]$s <- current.sample 
      }
      r.plate(plate.df)
    }
    
    session$resetBrush("plate_brush")
  })
  
  observeEvent(input$plate_dblclick, {
    req(input$plate_dblclick)
    df.row <- with(r.plate(), which(r == floor(input$plate_dblclick$y) & c == floor(input$plate_dblclick$x)))
    assign.sample(df.row, input$samp_select)
  })
  
  # Save the plate sample layout
  # https://stackoverflow.com/a/52909678/1193577
  observeEvent(input$save_layout, {
    showModal(modalDialog(
      tagList(
        textInput("save_layout_name", label = "Layout Name", placeholder = "")
      ), 
      title="Save sample layout",
      footer = tagList(actionButton("confirm_save_layout", "Save"),
                       modalButton("Cancel")
      )
    ))
  })
  
  observeEvent(input$confirm_save_layout, {
    req(input$save_layout_name)
    layout.nm <- gsub(" ", "_", trimws(input$save_layout_name))
    #TODO check name clash w/ existing layout
    
    if (!dir.exists("layouts")) dir.create("layouts", mode="0775")
    write.table(r.plate()[,c("r", "c", "w", "s")] %>% arrange(r, c),
                file.path("layouts", paste0(layout.nm, ".plt")),
                row.names=FALSE, quote=FALSE, sep="\t")
    
    removeModal()
  })
  
  # Load a saved plate sample layout
  observeEvent(input$load_layout, {
    showModal(modalDialog(
      tagList(
        selectInput("load_layout_name", label = "Layout Name", choices = gsub("\\.plt$", "", list.files("layouts/", "*.plt")))
      ), 
      title="Load sample layout",
      footer = tagList(actionButton("confirm_load_layout", "Load"),
                       modalButton("Cancel")
      )
    ))
  })
  
  observeEvent(input$confirm_load_layout, {
    req(input$load_layout_name)
    layout.path <- file.path("layouts", paste0(input$load_layout_name, ".plt"))
    new.layout <- read.table(layout.path, sep="\t", header=TRUE, colClasses = c("integer", "integer", "character", "character"))
    #TODO check no extraneous info in loaded file; check r, c columns too?
    
    # re-write s column; keep all other metadata associated with the wells
    if (all(r.plate()$w %in% new.layout$w)) {
      max.samp.n <- suppressWarnings(max(as.numeric(new.layout$s), na.rm=TRUE))
      updateNumericInput(session, "N_samps", value=max.samp.n)
      r.plate(
        left_join(new.layout, 
                  if ("s" %in% colnames(r.plate())) 
                    r.plate() %>% dplyr::select(-s)
                  else
                    r.plate()
                  ,
                  by=c("r", "c", "w"))
      )
    }
    
    removeModal()
  })
  
  # Display the plate
  output$plate <- renderPlot({
    p <- ggplot(r.plate()) +
      geom_rect(aes(xmin=c, xmax=c+1, ymin=r, ymax=r+1, fill=nchar(as.character(s))==0), color="black") +
      geom_text(aes(x=c, y=r, label=w, hjust="left"), nudge_x=0.12, nudge_y=-0.12) +
      geom_text(aes(x=c, y=r, label=s, hjust="right"), nudge_x=.9, nudge_y=-.7)
    
    if ("Assay" %in% colnames(r.plate())) {
      # file loaded
      p <- p + 
        geom_text(aes(x=c, y=r, label=short.name, hjust="left"), nudge_x=0.12, nudge_y=-0.3)
    }
    
    p <- p +
      scale_fill_manual(values=c(`FALSE`="#ffcccc", `TRUE`="white")) +
      scale_y_reverse() +
      theme_void() + 
      guides(fill=FALSE)
    
    p
  })
  
  # Display table of current data summary
  output$assay_table <- renderTable(r.assay.summary())
  
  # Display methy plot
  output$assay_plot <- renderPlotly({
    if (!is.null(r.assay.summary())) {
      to.plot <- r.assay.summary() %>% ungroup() %>% filter(!is.na(m.mean))
      if (nrow(to.plot) == 0) {
        ggplot() + geom_text(aes(x=0, y=50), label="No data to plot.") + theme_void()
      } else {
        ggplot(to.plot, aes(x=pos, y=m.mean, color=s, group=s)) +
          geom_point() +
          geom_line() +
          geom_errorbar(aes(ymin=m.mean-m.sd, ymax=m.mean+m.sd), width=0.2) +
          xlab("CpG") + ylab("Methylation") +
          ylim(c(0, 100)) +
          theme_bw() +
          theme(legend.title=element_blank())
      }
    }
  })
  
  get.excel <- function() {
    library(openxlsx)
    wb <- createWorkbook()
    addWorksheet(wb, "pyro_summary")
    
    summ <- get.pyro.summary()
    N.samps <- length(unique(summ$s))
    
    block.mean <- dcast(summ, Assay+pos~s, value.var="m.mean", fill=NA)
    block.sd   <- dcast(summ, Assay+pos~s, value.var="m.sd",   fill=NA)
    block.N    <- dcast(summ, Assay+pos~s, value.var="m.N",    fill=0)
    
    writeData(wb, 1, block.mean, startCol=1, startRow=2, colNames=TRUE)
    writeData(wb, 1, block.sd[,-1:-2], startCol=2+N.samps+1, startRow=2, colNames=TRUE)
    writeData(wb, 1, block.N[,-1:-2],  startCol=2+2*N.samps+1, startRow=2, colNames=TRUE)
    
    mergeCells(wb, 1, 3:(2+N.samps), 1)
    mergeCells(wb, 1, (2+N.samps+1):(2+2*N.samps), 1)
    mergeCells(wb, 1, (2+2*N.samps+1):(2+3*N.samps), 1)
    
    writeData(wb, 1, "Mean methylation", xy=c(3, 1))
    writeData(wb, 1, "Methylation SD",   xy=c(2+N.samps+1, 1))
    writeData(wb, 1, "Methylation N",    xy=c(2+2*N.samps+1, 1))
    
    centered.style <- createStyle(halign="center")
    addStyle(wb, 1, centered.style, rows=1, cols=c(3, 2+N.samps+1, 2+2*N.samps+1))
    bold.style <- createStyle(textDecoration = "bold")
    addStyle(wb, 1, bold.style, rows=c(1,2), cols=1:(2+3*N.samps), gridExpand=TRUE, stack=TRUE)
    addStyle(wb, 1, bold.style, rows=1:(nrow(block.mean)+2), cols=c(1,2), gridExpand=TRUE, stack=TRUE)
    setColWidths(wb, 1, 1, widths="auto")
    
    wb
  }
  
  # Save / download processed methy data as tsv
  output$export_tsv <- downloadHandler(
    filename = function() {
      paste0(gsub("(.*)\\.(.*)", "\\1", basename(r.input.file())), ".pyro_summary.tsv")
    },
    content = function(con) {
      write.table(get.pyro.summary(), con, row.names = FALSE, quote=FALSE, sep="\t")
    }
  )
  
  # Save / download processed methy data as tsv
  output$export_xlsx <- downloadHandler(
    filename = function() {
      paste0(gsub("(.*)\\.(.*)", "\\1", basename(r.input.file())), ".pyro_summary.xlsx")
    },
    content = function(con) {
      xlsx.wb <- get.excel()
      saveWorkbook(xlsx.wb, con)
    }
  )
  
  # Display info about well on hover
  # https://gitlab.com/snippets/16220
  output$hover_info <- renderUI({
    req(input$plate_hover)
    
    hover <- input$plate_hover
    hover.idx <- with(r.plate(), which(r == floor(hover$y) & c == floor(hover$x)))
    if (length(hover.idx) == 0) return(NULL)
    entry <- r.plate()[hover.idx,]
    #cat(str(entry),"\n")
    out.html <- paste0("<b> Well: </b>", entry$w, "<br/>",
                       "<b> Sample: </b>", entry$s, "<br/>")
    max.left <- 0.8
    
    if ("Assay" %in% colnames(entry)) {
      e.methy <- entry[1,endsWith(colnames(entry), "Methylation")]
      e.qual  <- as.character(unlist(entry[1,endsWith(colnames(entry), "Quality")]))
      valid.CG <- which(nchar(e.qual) > 1) # blank quality -> no CG in assay
      short.q <- substr(e.qual[valid.CG],1,1)
      bad.q <- which(short.q %in% c("F", "I"))
      if (length(bad.q) > 0) {
        short.q[bad.q] <- sprintf('<b><span style="color:red">%s</span></b>', short.q[bad.q])
      }
      pos.str <- paste(sprintf("  <b>%-3d</b>", valid.CG), collapse=" | ")
      methy.str <- paste(sprintf("%4.1f%%", unlist(e.methy[,valid.CG])), collapse=" | ")
      qual.str <- paste(sprintf("  %s  ", short.q), collapse=" | ")
      
      out.html <- paste0(out.html,
                         "<b> Assay: </b>", entry$Assay, "<br/>",
                         "<b> Methylation data: </b>", "<br/>",
                         "<pre>",
                         pos.str, "<br/>",
                         methy.str, "<br/>",
                         qual.str, "<br/>",
                         "</pre>")
      
      max.left <- 0.5
    }
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    #left_pct <- min(left_pct, max.left) # don't let it get too far to the right
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    side_str <- "left:"
    if (left_pct > 0.5) {
      side_str <- "right:"
      left_pct <- 1 - left_pct
    }
    side_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    #cat(left_pct, top_pct, left_px, top_px, "\n")
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.95); ",
                    side_str, side_px + 2, "px; top:", top_px + 2, "px;")
    
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(out.html))
    )
  })
  
}
