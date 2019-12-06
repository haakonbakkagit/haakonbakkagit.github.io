library(shiny)

ui<-fluidPage(
  shinyUI(navbarPage("My Application",
                     tabPanel("Tab 1",
                              sidebarLayout(
                                sidebarPanel(
                                  radioButtons("data_tab1", "Data:",
                                               c("Temperature" = 1,"Unemployment" = 2),
                                               selected="1",
                                               inline=TRUE)
                                ),
                                mainPanel("The graph of data:",
                                          plotOutput("plot_tab1"))
                              )
                              ),
                     tabPanel("Tab 2",
                              sidebarLayout(
                                sidebarPanel(
                                  radioButtons("data_tab2", "Data:",
                                               c("Temperature" = 1,"Unemployment" = 2),
                                               selected="1",
                                               inline=TRUE),
                                  
                                  conditionalPanel(
                                    condition = "input.data_tab2 == '1'",
                                    sliderInput(inputId = "theta1_tab2_1",
                                                label = "Choose a fixed value",
                                                value = 0, min = -3, max = 3)
                                  ),
                                  
                                  conditionalPanel(
                                    condition = "input.data_tab2 == '2'",
                                    sliderInput(inputId = "theta1_tab2_2",
                                                label = "Choose a fixed value",
                                                value = 0.5, min = -1.5, max = 2.5)
                                  )

                                ),
                                mainPanel("The fitting result corresponding different fixed values of θ1.",
                                          plotOutput("plot_tab2"))
                              )
                              ),
                     tabPanel("Tab 3", 
                              # titlePanel("The Fitting Result Corresponding Different Fixed Value"),
                              
                              sidebarLayout(
                                
                                sidebarPanel(
                                  sliderInput(inputId = "u_tab3",
                                              label = "Choose the value of u",
                                              value = 0.06, min = 0.0, max = 1.0),
                                  
                                  sliderInput(inputId = "a_tab3",
                                              label = "Choose the value of alpha ",
                                              value = 0.01, min = 0.0, max = 1.0)
                                ),
                                
                                mainPanel("The PC prior for precision is as following with different values of parameters (u,α).",
                                          plotOutput("plot_tab3"))
                              )
                              ),
                     tabPanel("Tab 4", 
                              # titlePanel("The Fitting Result Corresponding Different Fixed Value"),
                              
                              sidebarLayout(
                                
                                sidebarPanel(
                                  sliderInput(inputId = "u_tab4",
                                              label = "Choose the value of u",
                                              value = 0.9, min = 0.0, max = 1.0),
                                  
                                  sliderInput(inputId = "a_tab4",
                                              label = "Choose the value of alpha ",
                                              value = 0.9, min = 0.0, max = 1.0)
                                ),
                                
                                mainPanel("The PC prior for the correlation ρ with ρ = 1 as the base-model is as following with different values of parameters (u,α).",
                                          plotOutput("plot_tab4"))
                              )
                              ),
                     tabPanel("Tab 5", 
                              # titlePanel("The Fitting Result Corresponding Different Fixed Value"),
                              
                              sidebarLayout(
                                sidebarPanel(
                                  radioButtons("data_tab5", "Data:",
                                               c("Temperature" = 1,"Unemployment" = 2),
                                               selected="1",
                                               inline=TRUE),
                                  conditionalPanel(
                                    condition = "input.data_tab5 == '1'",
                                    sliderInput(inputId = "u_tab5_1",
                                               label = "Choose a value of u in pc.prec",
                                               value = 0.06, min = 0, max = 0.12)
                                  ),
                                  conditionalPanel(
                                    condition = "input.data_tab5 == '2'",
                                    sliderInput(inputId = "u_tab5_2",
                                                label = "Choose a value of u in pc.prec",
                                                value = 0.06, min = 0.03, max = 0.09)
                                  )
                                ),
                                mainPanel("The fitting result corresponding different fixed values of u in PC prior.",
                                          plotOutput("plot_tab5"))
                              )
                              )
)))

server<-function(input,output){
  
  output$plot_tab1<-renderPlot({
    if(input$data_tab1=="1"){
      
      load("../data/temperature-data")
      
      n <- length(tmed)
      data <-data.frame(y=tmed,t=1:n)
      
      pd <- pretty(c(dates, max(dates+30)), n=13)
      par(mfrow=c(1,1), mar=c(3,3,0.5,2), mgp=c(2,.7,0), las=2, xaxs='i')
      plot(dates, tmed, type='l', lwd=2,
           axes=FALSE, xlab='day', ylab='Temperature')
      abline(h=0)
      abline(h=3*(-8:9), v=pd, lty=3, col=gray(.5))
      box()
      axis(2, 3*(-8:9)); axis(4, 3*(-8:9))
      axis(1, pd, months(pd, TRUE))
      
    }else{
      
      temp = read.csv("../data/harmonised-unemployment-rates-mo.csv")
      n = nrow(temp)
      data = data.frame(y = temp[,2], t=1:n)
      dates <- temp[,1]
      
      plot(dates, data$y, lwd=2,
           xlab='month', ylab='Unemployment Rates')
      lines(dates,data$y)
      abline(h=2*(-8:9), lty=2, col=gray(.5))
      
    }
    
    
  })
  
  output$plot_tab2<-renderPlot({
    library(INLA); library(fields) 
    
    if(input$data_tab2=="1"){
      
      load("../data/temperature-data")
      
      n <- length(tmed)
      data <-data.frame(y=tmed,t=1:n)
      
      family <- "gaussian"
      
      hyper2 = list(theta1=list(initial=input$theta1_tab2_1, fixed=T))
      formula2 <- y~ f(t,model='ar1',hyper=hyper2)
      
      res2 <- inla(formula=formula2,data=data,family=family)
      plot(data$y, col="blue",
           ylab="fitting result")
      lines(res2$summary.random$t[ ,"mean"]+res2$summary.fixed$mean[1])
      
    }else{
      
      temp = read.csv("../data/harmonised-unemployment-rates-mo.csv")
      n = nrow(temp)
      data = data.frame(y = temp[,2], t=1:n)
      dates <- temp[,1]
      
      family <- "gaussian"
      
      hyper2 = list(theta1=list(initial=input$theta1_tab2_2, fixed=T))
      formula2 <- y~ f(t,model='ar1',hyper=hyper2)
      
      res2 <- inla(formula=formula2,data=data,family=family)
      plot(data$y, col="blue",
           ylab="fitting result")
      lines(res2$summary.random$t[ ,"mean"]+res2$summary.fixed$mean[1])
      
    }
    
  })
  
  output$plot_tab3<-renderPlot({
    lambda=-log(input$a_tab3)/input$u_tab3
    m=paste("lambda=",lambda)
    prec=(0:10000); plot(prec, main=m, ylab="Probability density", xlab= "prec" , type="l", inla.pc.dprec(prec, input$u_tab3, input$a_tab3, log = FALSE))
  })
  
  output$plot_tab4<-renderPlot({
    
    library("BB")
    fun <- function(x) { 
      f <- numeric(length(x)) 					
      f[1] <- exp(-x[1]*sqrt((1-input$u_tab4)))-input$a_tab4+input$a_tab4*exp(-sqrt(2)*x[1])
      f 
    } 
    startx <- c(0.3)
    result = dfsane(startx,fun,control=list(maxit=2500,trace = FALSE))
    lambda = result$par
    m=paste("lambda=",lambda)
    
    cor=(500:1000)/1000; plot(cor, main=m, ylab="Probability density", xlab= "cor" , type="l", inla.pc.dcor1(cor, input$u_tab4, input$a_tab4, log = FALSE))
  })
  
  output$plot_tab5<-renderPlot({
    library(INLA); library(fields) 
    
    if(input$data_tab5=="1"){
      
      load("../data/temperature-data")
      
      n <- length(tmed)
      data <-data.frame(y=tmed,t=1:n)
      
      family <- "gaussian"
      
      hyper4 <- list(theta1 = list(prior="pc.prec", param=c(input$u_tab5_1, 0.01)),
                     theta2 = list(prior="pc.cor1", param=c(0.9, 0.9)) )
      formula4 <- y~ f(t,model='ar1',hyper=hyper4)
      res4 <- inla(formula=formula4,data=data,family=family)
      
      plot(data$y, col="blue",ylab="fitting result")
      lines(res4$summary.random$t[ ,"mean"]+res4$summary.fixed$mean[1])
      
    }else{
      
      temp = read.csv("../data/harmonised-unemployment-rates-mo.csv")
      n = nrow(temp)
      data = data.frame(y = temp[,2], t=1:n)
      dates <- temp[,1]
      
      family <- "gaussian"
      
      hyper4 <- list(theta1 = list(prior="pc.prec", param=c(input$u_tab5_2, 0.008)),
                     theta2 = list(prior="pc.cor1", param=c(0.9, 0.9)) )
      formula4 <- y~ f(t,model='ar1',hyper=hyper4)
      res4 <- inla(formula=formula4,data=data,family=family)
      
      plot(data$y, col="blue",ylab="fitting result")
      lines(res4$summary.random$t[ ,"mean"]+res4$summary.fixed$mean[1])
      
    }
    
  })
  
}
shinyApp(ui = ui, server = server)








