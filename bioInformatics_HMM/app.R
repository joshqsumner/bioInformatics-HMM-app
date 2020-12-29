library(DT)
library(tidyverse)
library(plotly)
library(stringr)
library(RColorBrewer)
library(shiny)
library(scales)
library(shinydashboard)
library(lubridate)
library(grDevices)
library(HMM)
library(seqHMM)
library(depmixS4)
library(aphid)
library(ape)
library(shiny)

# outside of visible app
## HMM Tab's hmm
States = c("Fair", "Unfair") #two dice
Symbols = 1:6 #both are six sided dice
transProbs = matrix(c(0.99, 0.01, 0.02, 0.98), c(length(States), 
                                                 length(States)), byrow = TRUE) #transition probabilities as a 2x2 matrix (2 states)
emissionProbs = matrix(c(rep(1/6, 6), c(rep(0.1, 5), 0.5)),
                       c(length(States), length(Symbols)), byrow = TRUE) #emission probs as 6x6 matrix (6 sides)
HMM_hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs)

## aphid tab's hmm
accessionNumbers <- c("U15717", "U15718", "U15719", "U15720",
                      "U15721", "U15722", "U15723", "U15724") #some tanager birds, who doesn't love birds???
birds <- ape::read.GenBank(accessionNumbers) #if you don't love birds I don't love you
#jk you're probably alright either way. they're my annotations, I'll do what I want.
#okay I should mention more about this, if you get an error running this app and it's something about internet/a website
#being unreachable then it is the read.GenBank function. It needs to be able to reach ncbi and download things.
#that can be finnicky if you are like me and have bad wifi. If it fails you can try the URL in a browser to see what happens but
#it'll probably fail too. NCBI will sometimes be down (https://www.ncbi.nlm.nih.gov/down/) so that will also throw a wrench in this
#i'll try to make it fail gracefully if possible, but it might be beyond the scope of what i want to do for this project.
birds3<-birds[1:3]
birdsR<-birds[4:8]
birds.PHMM <- aphid::derivePHMM(birds3) #if you thought they had to be helpful that's just too bad.

#this is all for the seqHMM tab, and that is part of why I am not a fan
data("mvad", package = "TraMineR") #mvad data from traminer, which is another package used for sequencing
mvad_alphabet <- c("employment", "FE", "HE", "joblessness", "school", "training")
mvad_labels <- c("employment", "further education", "higher education",
                 "joblessness", "school", "training")
mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad_seq <- seqdef(mvad, 15:86, alphabet = mvad_alphabet, states = mvad_scodes, labels = mvad_labels, xtstep = 6)
# seqdef is from TraMineR
# basically it makes a special kind of sequence object called an stslist.
attr(mvad_seq, "cpal") <- colorpalette[[123]] #who doesn't enjoy color palette number 123? Nobody, that's who.

# emission probs start point
emiss <- matrix(
  c(0.05, 0.05, 0.05, 0.05, 0.75, 0.05, # SC
    0.05, 0.75, 0.05, 0.05, 0.05, 0.05, # FE
    0.05, 0.05, 0.05, 0.4,  0.05, 0.4,  # JL, TR
    0.05, 0.05, 0.75, 0.05, 0.05, 0.05, # HE
    0.75, 0.05, 0.05, 0.05, 0.05, 0.05),# EM
  nrow = 5, ncol = 6, byrow = TRUE)

# first values for transition probs
trans <- matrix(0.025, 5, 5)
diag(trans) <- 0.9

# first values for initial probabilities
initial_probs <- c(0.2, 0.2, 0.2, 0.2, 0.2)

init_hmm_mvad <- build_hmm(observations = mvad_seq,
                           transition_probs = trans, emission_probs = emiss,
                           initial_probs = initial_probs)

set.seed(21)
fit_hmm_mvad <- fit_model(init_hmm_mvad, control_em = list(restart = list(times = 2)))
hmm_mvad <- fit_hmm_mvad$model #The fit_hmm_mvad object has a few pieces

# Okay this is also seqHMM but it's going to be a simpler version to explain the graph better
rd<-mvad #make a new dataframe from the mvad data to mess with
varList<-colnames(rd)[15:86] #take the relevant columns
for (i in varList) { #iterate through those columns
  rd[[i]] <- factor(ifelse(rd[[i]] %in% c("FE", "school", "training", "HE"), "education", as.character(rd[[i]]))) #simplify them
}
rd_alphabet <- c("education", "employment", "joblessness")
rd_labels <- c("education", "employment", "joblessness")
rd_scodes <- c("ED", "EM", "JL")
rd_seq <- seqdef(rd, 15:86, alphabet = rd_alphabet, states = rd_scodes, labels = rd_labels, xtstep = 6)
attr(rd_seq, "cpal") <- colorpalette[[5]]
rdemiss <- matrix(
  c(0.15, 0.3, 0.55, 
    0.5, 0.4, 0.1, 
    0.65, 0.25, 0.1),
  nrow = 3, ncol = 3, byrow = TRUE)

rdtrans <- matrix(0.05, 3, 3)
diag(rdtrans) <- 0.9
initial_probs <- c(1/3, 1/3, 1/3)
init_hmm_rd <- build_hmm(observations = rd_seq,
                         transition_probs = rdtrans, 
                         emission_probs = rdemiss,
                         initial_probs = initial_probs)

set.seed(21)
fit_hmm_rd <- fit_model(init_hmm_rd, control_em = list(restart = list(times = 1)))
#end of seqHMM stuff

# depmixS4 outside app stuff
data(speed)
data(sp500)


# Define UI
ui <- dashboardPage(skin="black",
  dashboardHeader(title = "BioInformatics M21-550"),
  dashboardSidebar(sidebarMenu(
    menuItem("Home",
             tabName = "home",
             icon = icon("readme")),
    menuItem("HMM::", tabName = "hmm_tab", icon=icon("chart-pie")), #these tabs will be for software/package options
    menuItem("aphid::", tabName = "aphid_tab", icon=icon("clipboard-list")),
    menuItem("seqHMM::", tabName = "seqhmm_tab", icon=icon("clipboard-list")), 
    menuItem("depmixS4::", tabName = "depmix_tab", icon=icon("award")),
    menuItem("vanillaICE::", tabName = "VI_tab", icon=icon("book")),
    menuItem("HMMER (CLI)", tabName = "hmmer_tab", icon=icon("cogs")), 
    menuItem("Other", tabName = "other_tab", icon=icon("globe-americas")), 
    menuItem("Other2", tabName = "other_tab_2", icon=icon("address-card")))), #I need to go update these icons, these ones are just random
  dashboardBody(
    tabItems(
      # Home Tab
      tabItem("home", 
              tabBox(
                id = "tabset1", height = "250px", width="600px",
                tabPanel("Overview", tags$ol(
                  tags$li("The tabs to the left show examples of different HMM softwares/R packages."),
                  tags$li("Code can be downloaded here: ."), 
                  tags$li("Presentation can be downloaded here: ")
                )),
                tabPanel("Finding code/resources", tags$ol(
                  tags$li("Code can be downloaded here: ."), 
                  tags$li("Data can be viewed in tabs and is annotated in github"), 
                  tags$li("Presentation: ")
                )))),
      #R packages Tabs
      
      tabItem("hmm_tab", 
              tabBox(title="R HMM Package", height="1000px", width="1000px",
                     tabPanel("How HMM:: works", width="800px", p("For get started with HMM::
                               this package provides a toy example through the", 
                                                                  code('dishonestcasino()'), "function.
                               It takes no arguments and will run in the R console where 
                               you can just click through to show this example model."),
                       img(src="DCex.png"),
                       br(),
                       br(),
                     p("\nThat will output a base R plot similar to the one on the next tab."),
                     br(),
                     img(src="DCex2.png"),
                     br(),
                     p("\n The next tab is using different code to make a similar visualization. \n"),
                     br(),
                     p("The key function starting off is ", code('initHMM()'), "which makes a 
                        discrete hidden markov model from a vector of states and emissions 
                       (symbols) as well as an initial probability vector, transition probability matrix,
                        and emission probability matrix. That will return an HMM object, which has 5 elements in it:"),
                     tags$ol(
                       tags$li("States Vector with the names of the states."),
                       tags$li("Symbols Vector with the names of the symbols (emitted symbols)."),
                       tags$li("startProbs Annotated vector with the starting probabilities of the states."),
                       tags$li("transProbs Annotated matrix containing the transition probabilities between the states."),
                       tags$li("emissionProbs Annotated matrix containing the emission probabilities of the states.")
                       ),
                    p(" With an initialized HMM there are a variety of
                       functions to find ", code('forward()'), ", ", code('posterior()'), ", and ", code('backward()'), " functions compute the respective probabilities.
                      As well as ", code('simHMM()'), " and ", code('viterbi()'), " functions which allow the use/simulation with the initialized HMM.
                      HMM:: also has a ", code('baumWelch()'), " function which finds locally optimized parameters for the HMM given a set of data.
                      The same task can also be done with ", code('viterbiTraining()'), " which can basically be thought of as fast and loose ", code('baumWelch()'))
                     
                     ),
                     tabPanel( 
                       title="Dishonest Casino HMM",
                       numericInput('die_throws', 'Number of Rolls',
                                   value = 1000,
                                   min = 2,
                                   max = 1000000),
                       plotOutput("dishonestCasinoBase", width="750px", height="500px"),
                       br(),
                       br(),
                       br(),
                       p("This plot is being made by initializing a static hmm with the probabilities we used in class 
                         (code snippet below for reference), simulating N paths using ", code('simHMM(hmm, N)'), " then finding the viterbi path using
                         ", code('viterbi(hmm, simulatedHMM$observation)'), "."),br(),
                       code('Dice = c("Fair", "Unfair")'), br(),
                       code('Outcomes = 1:6'),br(),
                       code('tprobs = matrix(c(0.99, 0.01, 0.02, 0.98), c(length(Dice), length(Dice)), byrow = TRUE)'),br(),
                       code('eprobs = matrix(c(rep(1/6, 6), c(rep(0.1, 5), 0.5)), c(length(Dice), length(Outcomes)), byrow = TRUE)'),br(),
                       code('hmm = initHMM(Dice, Outcomes, transProbs = tprobs, emissionProbs = eprobs)')
                       ), 
                     tabPanel( 
                       title="Dishonest Casino Resample",
                       numericInput('die_throws_2', 'Number of Rolls per Trial',
                                    value = 100,
                                    min = 2,
                                    max = 10000),
                       numericInput('casino_trials', 'Number of Trials (max 500)',
                                    value = 100,
                                    min = 1,
                                    max = 500),
                       checkboxInput(
                         inputId = 'whichPlotter', label = 'Use Plotly', value = F
                       ),
                       conditionalPanel(
                         condition = 'input.whichPlotter == true', #for classmates, true is lowercase bc this is technically javascript, my code annotations will sometimes be helpful!
                         plotlyOutput('dishonestCasinoBootPlotly')
                       ),
                       
                       conditionalPanel(
                         condition = 'input.whichPlotter == false', #second verse same as the first
                         plotOutput('dishonestCasinoBoot', width="750px", height="500px")
                       ),
                       br(),
                       p("The HMM tends to do a very good job with this kind of task. Generally the averages are in the 90%+ range.
                         When you simulate smaller samples though (100 trials of 10-ish die rolls at a time, for instance) you can get a lot worse accuracy.
                         That is because it is very costly (from a probability standpoint) to change hidden states and the initial probabilities are 
                         very unbalanced (in the spirit of the example it is unlikely a casino starts off by robbing you, hopefully you'd catch on quickly).
                         In those situations a HMM will not work as well since the penalty for being wrong about the hidden state
                         is not outweighed on any path by the likelihood of the observed outcomes, even if they are all wrong. Increase your rounds playing baccarat during each trial,
                         or whatever you prefer to lose money at, and the average accuracy goes up quickly. With ~100 dice throws per trial it is very rare
                         to see any observations below ~45% accurate.")
                     ))),
      
      tabItem("aphid_tab", 
              tabBox(title = "R aphid Package",
                     id = "aphidTab", height="1000px", width="1000px",
                     tabPanel("How aphid:: Works",
                              p("This class were briefly introduced to aphid:: in lab 4. We mostly read documentation and created a
                                visual representation of an HMM like this:"),
                              plotOutput("aphidBasePlot"),
                              br(),
                              p("Aphid:: supports both HMMs and PHMMs and unlike ", code('HMM::'), " it is made with biological applications in mind and is essentially
                                based off of/made to interact with 'HMMER' software, which was developed at WashU by the author of the review paper I read to start this project. 
                              The package still has test data included for the classic dishonest casino example but also includes example globin data
                                (", code("globins"),") which allows for more relevant examples and learning. Before getting to that though, there 
                                are a couple of ways to initiate an HMM/PHMM in aphid:"),
                              tags$li(code("deriveHMM()"), " is a function to build an HMM from a list of sequences.
                                      Those sequences should be character vectors of emitted observations and the vectors 
                                      should be named to represent the hidden state from which each residue was emitted. 
                                      Importantly you can use DNAbin or AAbin classes from ", code("ape::"), " which saves effort 
                                      and potential confusion in converting those sequences into long character vectors. Any of those inputs will yield an HMM class object.", br(), 
                                      "The only required argument is the list of sequences, x, but there are more optional parameters", 
                                      code('deriveHMM(x, seqweights = NULL, residues = NULL, states = NULL, modelend = FALSE, pseudocounts = "background", logspace = TRUE)')),
                              tags$li(code("derivePHMM()"), " this function is the complement of ", code("deriveHMM()"), " which makes a profile HMM object of class 'PHMM'.
                                      ", code("derivePHMM()"), " has a large number of optional arguments but the most immediate difference is that the 'x' argument
                                      requires a matrix of aligned sequences or a list of unaligned sequences. These can be characters or the same special ",
                                      code("ape::"), " classes. Other interesting arguments include 'seqweight' to specify different weights for the sequences,
                                      'threshold' which lets you limit how much of the alignment can be gap space, and 'refine' which lets you change 
                                      between methods for refining parameters (such as viterbi)."),
                              tags$li(code("readPHMM()"), " is a helpful function especially if you have some experience with 'HMMER' or are working
                                      with a team who does. This function lets you read in a PHMM object from a different analysis and store it in your R session
                                      as a PHMM object. There is a companion ", code("writePHMM()"), " function which does the opposite and saves a PHMM object as text
                                      readable by 'HMMER' or ", code("readPHMM()"), " if you are so inclined."),
                              br(),
                              p("With a HMM/PHMM object ready you can plot it to get a general feel for what is going on or use 
                                other functions for predictions or evaluations which I'll go more into here?")
                              
                              ),
                     tabPanel("Example Data", 
                              p("The ", code("ape::"), " package has a set of 15 woodmouse 
                                DNA sequences which, using ", code("derivePHMM()"), " makes the following Profile HMM. Here all the possible options
                                for a DNA sequence are present (there are only the 4 nucleotide options, so it'd be really strange to not have them all shown)
                                so you do not need to specify a 'residue' argument."),
                              br(),
                              plotOutput("woodmousePlot"),
                              br(),
                              p("This visualization is coming from woodmouse DNA sequences in a DNAbin object in R.
                                Below is a printout of the first few positions of that sequence as a Table."),
                              br(),
                              dataTableOutput("woodmouseData"),
                              br(),
                              p("There is also a set of sample protein data in the ", code("aphid::"), " package."),
                              br(),
                              plotOutput("globinPlot"),
                              br(),
                              dataTableOutput("globinData"),
                              br(),
                              p("You may notice there are more emission options shown than are present in this little bit of sample data.
                                That is because the ", code("derivePHMM()"), " function's 'residues' argument tells the model what the possible emission states
are depending on what the sequence is coming from. The default is NULL which only reads emission states from the input data
                                , but you can use a character vector of 'levels' or a set string from 'RNA', 'DNA', 'AA', or 'AMINO'."),
                              br(),
                              p("These are nice graphics but ", code("aphid::"), " has ways to directly use the PHMM/HMM you've made as well. 
                                The ", code("align()"), "function can take a list of unaligned sequences (as DNAbin/AAbin objects or a character vectors and a model to 
                                align the sequences with. Here is an example of how that works with the short globin data and woodmouse data."),
                              br(),
                              code("data(globins)"),
                              br(),
                              code("globins.PHMM <- derivePHMM(globins, residues = 'AMINO'"),
                              br(),
                              p(code("globinsUnaligned <- unalign(globins)"), "This step is only necessary because I'm using the sample data which comes aligned"),
                              br(),
                              code("globinsReAligned<-align(globins2, model = globins.PHMM, seqweights = NULL, residues = 'AMINO')"),
                              br(),
                              p("Now ", code('globinsUnaligned'), " is a list which contains:"),
                              code("$HBA_HUMAN"),
                              br(),
                              code('[1] "V" "G" "A" "H" "A" "G" "E" "Y"'),
                              br(),
                              code("$HBA_HUMAN"),
                              br(),
                              code('[1] "V" "N" "V" "D" "E" "V"'),
                              br(),
                              code("$MYG_PHYCA"),
                              br(),
                              code('[1] "V" "E" "A" "D" "V" "A" "G" "H"'),
                              br(),
                              code("$GLB3_CHITP"),
                              br(),
                              code('[1] "V" "K" "G" "D"'),
                              br(),
                              code("$GLB5_PETMA"),
                              br(),
                              code('[1] "V" "Y" "S" "T" "Y" "E" "T" "S"'),
                              br(),
                              code("$LGB2_LUPLU"),
                              br(),
                              code('[1] "F" "N" "A" "N" "I" "P" "K" "H"'),
                              br(),
                              code("$GLB1_GLYDI"),
                              br(),
                              code('[1] "I" "A" "G" "A" "D" "N" "G" "A" "G" "V"'),
                              br(),
                              p("And the ", code("globinsReAligned"), " object contains this aligned sequence as a matrix/array:"),
                              br(),
                              dataTableOutput("globinsReAligned")
                     ),
                       tabPanel("GenBank Data", 
                                p("Here is an exmaple of making and using an ", code("aphid::"), " PHMM with other data. For this example
                                  I am also using the ", code("ape::"), " package that the woodmouse data came from. The ", code("ape::"),
                                  " package has a lot of great tools for doing this sort of thing and I highly recommend exploring it further but I'm 
                                  not going to get into many details here besides that I'm using the ", code("read.GenBank()"), "with a list of accession numbers
                                  to pull data from GenBank for this example and automatically make it into a DNAbin object. Full disclosure, I'm using accession numbers from a vignette
                                  since I do not really have another group of sequences in mind to explore. An important note about this is that NCBI will sometimes
                                  be down for maintenance or your wifi could be problematic and either of those will keep this from working and are just out
                                  of your control. That was a problem multiple times while finishing this so it is possible there will be errors for you reproducing it 
                                  just based on timing. I have no idea how consistently it works in China but I am curious if you end up trying."),
                                br(),
                                p(code('accessionNumbers <- c("U15717", "U15718", "U15719", "U15720", "U15721", "U15722", "U15723", "U15724")'), br(),
                                  code('birds <- read.GenBank(accessionNumbers)'), br(),
                                  code('cbind(attr(birds, "species"), names(birds))'), br(),
                                  "This yields a DNAbin object with these 8 sequences in it. That's great, but to make a model we want a reference so 
                                  I'm going to take the first three sequences out and use that to make a model to align the others.", br(),
                                  code("birds3<-birds[1:3]"), br(),
                                  code('birds.PHMM <- aphid::derivePHMM(birds3)'),
                                  br(), 
                                  code('plot(birds.PHMM, from = 0, to = 10, main = "Three Birds profile HMM")'), br(),
                                  "Now we have a working PHMM for this data. These sequences are >1000 bases long, so I'm just going to print out the first
                                  ten positions of the PHMM but you get the point."),
                                plotOutput("firstThreeBirdsPHMM"),
                                p("The model can then be used to align sequences.", br(), 
                                  code("birdsR<-birds[4:8]"), br(),
                                  code("align(birdsR, model=birds.PHMM)"), br(),
                                  "Here I just used the first three bird's DNA to make the model, then aligned the remaining 5 birds using the model."),
                                textOutput("alignedBirds"),
                                br(),
                                p("Sorry Ryan, I was going to add an option to input your own accession numbers but I decided I'd rather not.")
                                )
                  
                     )),
      
      
      tabItem("seqhmm_tab", 
              tabBox(title="seqHMM:: R Package", height="1000px", width="1000px",
                     tabPanel(
                       title="How it works",
                       p(code("seqHMM::"), " contains a variety of functions for building a variety of markov models including
                         hidden markov models, mixture hidden markov models, markov models, mixture markov models, and latent class models. 
                          Each kind of model can be initiated with a ", code("build_*()"), 
                         " function (", code("build_hmm()"), " for example). These functions take somewhat different arguments due to their
                         respective roles and parameters, I only plan to work with the ", code("build_hmm()"), " function here."), br(),
                       p("The ", code("build_hmm()"), " function takes observations and n_states arguments and has a variety of optional arguments for 
                         setting transition, initial, and emission probabilities. 
                         The next tab will have an example of using this function along with some help from another package
                         (", code("TraMineR::"),") to generate an HMM.")
                     ),
                     tabPanel("Example HMM", 
                              p("This HMM is going to use mvad data from ", code("TraMineR"), " as well as some helpful functions
                                from the same package. Just like how ", code("aphid"), " and ", code("ape"), "work well together,
                                ", code("seqHMM"), " and ", code("TraMineR"), " are made to be used together. The mvad data is from a study of how people
                                transition from school to work, it looks like this but with 87 columns:"),
                              dataTableOutput("mvadHead"),
                              br(),
                              p("So we basically can see what is going on now, there are a bunch of columns of 'states' from each time
                                measurement as the study follows the subjects through their professional development. The levels of those states
                                are ", code("employment, further education, higher education, joblessness, school,"), " and ", code("training."), 
                                "The first step is to explicitly define those levels for our model in whatever format we want to use: ", br(),
                                code('mvad_alphabet <- c("employment", "FE", "HE", "joblessness", "school", "training")'), br(),
                                code('mvad_labels <- c("employment", "further education", "higher education", "joblessness", "school", "training")'),
                                br(),
                                code('mvad_scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")'), br(),
                                "So now we have a few lists of states and labels and we need to combine those with our data to make a special ",
                                code("stslist"), "object since that is the format required for ", code("seqHMM::build_hmm()"), ". To do that we will use ",
                                code("TraMineR::seqdef()"), " where we will need to say what data to use, which columns contain relevant sequences, what alphabet, states, and labes to use
                                , and what kind of time gaps to use (the xtstep argument):",
                                code("mvad_seq <- seqdef(mvad, 17:86, alphabet = mvad_alphabet, states = mvad_scodes, labels = mvad_labels, xtstep = 6)"), br(),
                                "Now ", code("mvad_seq"), " contains an stslist which is a list with 14 elements including a dataframe. The stslist is actually
                                a pretty interesting object but for now I'm just going to mention that there are methods for using basic functions on it or use TraMineR specific versions:",
                                code("plot(mvad_seq)")),
                              plotOutput("mvad_seq_plot"),
                              br(),
                              p("That's one of the situations where having more control over the plotting would be nice, but it's beside the point for now.
                                Next we need to define the HMM with ", code("build_hmm()"), " passing the stslist and matrices of emission/transition probabilities
                                and a vector of initial probabilities:", br(),
                                code("init_hmm_mvad <- build_hmm(observations = mvad_seq, transition_probs = trans, emission_probs = emiss, initial_probs = initial_probs)"), br(),
                                "That yields an hmm object with whatever probabilities we defined in those matrices/vectors. Now we can use the ", code("seqHMM::fit_model()"), 
                                " function to estimate the parameters and return a list with the model and different results from it.", br(),
                                code("fit_hmm_mvad <- fit_model(init_hmm_mvad, control_em = list(restart = list(times = 2)))"),
                                br(),
                                "The 'times' argument is telling the expectation maximization algorithm how many times to restart from initial random values. Normally you would use
                                more than 2 here, but for the sake of time and since we aren't doing anything too important with this 2 is fine. You can also specify
                                in that list of arguments whether you want the transition/emission probabilities to vary on each iteration among other things. There is a lot of flexibility with
                              ", code("fit_model()"), " if you have something specific you want it to do. Now we can display the model, although seqHMM does not take great pains
                                to make the visualizations easily interpretable/use good visualization practices like the other packages we have checked so far do (although they do bother to provide a list of 200 color palettes
                                which they think are good for showing these models which I think is an odd but fun choice, I arbitrarily choose #123):", br(), 
                                code("plot(fit_hmm_mvad$model)")),
                              plotOutput("mavdHMMPlot"), br(),
                              p("I really want to like this as a plot, but it's pretty difficult to intuitively read it, 
                                pie charts are a terrible idea just in general, the text is not repelled nor is it trivial to repel it, the hmm object is not particularly conducive to ggplot/plotly, etc.
                                Each pie corresponds to a hidden state with the emission states probabilities shown as slices of the pie. Transition probabilities
                                are shown with arrows between different pies. You get a ton of versatility from this package,
                                but most of the time it really seems like ", code("HMM"), " or ", code("aphid"), " would work just as well while going a lot further to make your code
                                easily readable and provide you with a straightforward method to get some visualizations to explain the HMM should that be important. In practical settings for sequence
                                analysis, visualization is the least of your concerns though so ", code("seqHMM"), " is still a very useful package for bioinformatics, it just might not lend itself to 
                                examples as well. To try to show what is happening a little easier I made a even simpler toy example:"),
                              plotOutput("rdHMMPlot"),
                              br(),
                              p("There are equal initial probabilities, these are the emission probabilities:"),
                              br(),
                              dataTableOutput("rdEmiss"),
                              p("And these are the transition probabilities:"),
                              br(),
                              dataTableOutput("rdTrans")
                              
                              )#,
                    # tabPanel("Thoughts on seqHMM", "If i have other thoughts maybe those go here?")
                    # )
              )),  
      
      tabItem("depmix_tab", 
              tabBox(title="depmixS4:: R Package", height="1000px", width="1000px",
                     tabPanel(
                       title="How depmixS4 works",
                      p("The ", code("depmixS4"), " package is a relatively new R package for handling 'dependent mixture models' (which is debatably the same as an HMM).
                      A mixture model is just a probability model for measuring subpopulations within a greater population which does not require that the subpopulations are already
                        labelled in your data. For practical purposes it is generally fine to consider HMMs as a kind of mixture model where the distributions of emission probabilities are different between 
                        hidden states and the total model could be considered a mixture model of multiple distributions. Math aside, this package uses 'dependent mixture model' and HMM interchangably."), br(),
                      p(code("depmixS4::depmix()"), "is one of the primary functions of interest and is used to create a special ",
                        code("depmix"), " object. That object has a variety of 'slots' with information describing the model but it is not a completed model yet and the parameters are unfitted.
                        The model can be fit to your data with the ", code("depmixS4::fit()"), " package. The ", 
                        code("depmix()", " function provides a lot of versatility to how you want to build your model, the next tab has 
                             several examples of how it can be used in different applications."))
                     ),
                     tabPanel("Speed Data", 
                              p("Using a sample dataset of response times given the presence or absence of incentivization we can make an HMM.
                                Pick a number or results to resample (original data is fairly small, resampling will work better to show HMM properties),
                                a number to train on (remainder will be shown in test set), and set a seed if you want to."),
                              
                              numericInput("depTraining", "Set Training Data Size", value=300, min=10, max=1000),
                              numericInput("depTest", "Set Test Data Size", value=600, min=10, max=10000),
                              checkboxInput("depSeedCheck", "Set a Sampling Seed?", value=T),
                              conditionalPanel(condition = "input.depSeedCheck == 'true'",
                                               numericInput("depSeed", "Specify Seed:", value=123, min=1, max=10000)
                              ),
                              checkboxInput(
                                inputId = 'depWhichPlotter', label = 'Use Plotly?', value = F
                              ),
                              conditionalPanel(
                                condition = 'input.depWhichPlotter == true', 
                                plotlyOutput('depPlotly')
                              ),
                              conditionalPanel(
                                condition = 'input.depWhichPlotter == false',
                                plotOutput('depGGplot', width="750px", height="500px")
                              ),
                              p("Playing with these results you should see some patterns that having larger data allows for better results. That is normal, but here it gets
                                another nuance which is that the amount of training data and test data are both equally influential (viterbi works well with lots
                                of test observations, poorly with smaller data) just like in the previous pages. You'll also probably notice that the results aren't that great.
                                That is a function of the data not having especially strong signs for which hidden state an emission comes from, considering the setting that makes some sense.")),
                      tabPanel("Stock Data",
                              p("Next will be some stock data from the ", code("sp500"), " dataset. Here is a look at the data:"),
                              br(),
                              dataTableOutput("sp500DT"),
                              br(),
                              p("To get a better feel for it, here are some time series plots of this data."),
                              plotOutput("sp500TSClose"),
                              br(),
                              p("And one of the log returns."),
                              br(),
                              plotOutput("sp500TSLR"),
                              p("This actually goes to show one of the things I like about ", code("depmixS4"), " which is that it works great
                                with very standard objects in R. Dataframes are fine, tidy manipulations are fine, no special object class needed for input. 
                                It might be a little lazy for me to have that as a high priority, but I like what I like and not having to look up how to clean/maneuver some 
                                weird bioconductor object saves me time and energy. So now I'll do something with the stock data."),
                              br(),
                              p("One interesting thing to do with this is to look at the log returns and try to categorize different market stages based on that.
                                It's not what we normally do but it's a cool example and actually the first place I learned anything about hidden markov models.",
                                code("depmixS4::depmix()"), " makes this really easy.", br(), 
                                code("mod_sp <- depmix(logret~1, nstates=2, data=sp500)"), br(),
                                code("fitted_mod_sp <- fit(mod_sp)"), br(),
                                code('plot(ts(depmixS4::posterior(fitted_mod_sp)[,2], start=c(1950,1), deltat=1/12), ylab="probability",
                                    main="Posterior probability of state 1 (volatile, negative markets).", frame=FALSE)'),br(),
                                "That code is first making a depmix object out of a dataframe, a formula, and a number of hidden states.
                                Next it is fitting the model which optimizes the (randomly initialized) parameters. Last it is converting
                                the posterior probabilities (which are returned as a dataframe) for state 1 into a time series starting from 1950 and plotting the probability of a volatile
                                market at any particular time based on the log returns. Comparing the probabilities from the model to the graph of 
                                log returns you can see the alignment to stock market crashes/booms (1987, for example). If you have multiple models, for instance
                                you are trying out multiple options for how many hidden states there should be to best fit some data, you can compare results
                                with the ", code("stats::AIC()"), "/", code("stats::BIC()"), " functions. AIC/BIC is also available in the ", 
                                code("summary()"), " output for a depmix object. You can get a bit of a feel for how that works below with some mock poisson data
                                which should show two hidden states being the best model, but you can confuse it with similar enough distributions."),
                              br(),
                              numericInput("lambda1", "Lambda for first part of data", value=1, min=0, max=5),
                              numericInput("length1", "Length of first part of data", value=50, min=1, max=5000),
                              numericInput("lambda2", "Lambda for second part of data", value=2, min=0, max=5),
                              numericInput("length2", "Length of second part of data", value=50, min=1, max=5000),
                              checkboxInput("bicSeedCheck", "Set a Sampling Seed?", value=T),
                              conditionalPanel(condition = "input.bicSeedCheck == 'true'",
                                               numericInput("depSeed", "Specify Seed:", value=456, min=1, max=10000)
                              ),
                              br(),
                              plotOutput("poissonBIC"),
                              p("Those points are being generated by checking the BIC of different models made like this for X in 1 through 5:", 
                                code("mX <- depmix(y~1, nstates=X,family=poisson(), data=ydf)", br(),
                                      "set.seed(1)", br(),
                                     "fmX <- fit(mX,em=em.control(maxit=200))"))
                              ),
                     tabPanel("Biological Data", 
                              p(""))
                     )),
      tabItem("VI_tab", 
              tabBox(title="vanillaICE:: BioConductor R Package", height="1000px", width="1000px",
                     tabPanel("Thoughts",
                       p("Not a fan.")
                       )
                     )),
      
      # HMMER Tab
      
      tabItem("hmmer_tab", 
              tabBox(title = "HMMER Software",
                     id = "tabset1", height = "250px",
                     tabPanel("History",  tags$ol(
                       tags$li("HMMER history"), 
                       tags$li("Washu Role"), 
                       tags$li("Sean Eddy")),
                     tabPanel("Use", tags$ol(
                       tags$li("Using HMMER"), 
                       tags$li("Steps etc"), 
                       tags$li("More steps/data needs")))))),
      
      ## Viterbi Algorithm should be included since that comes up a bunch
      tabItem("Viterbi Paths", 
              tabBox(title = "other",
                     id = "tabset1", height = "250px", width="600px",
                     tabPanel("Viterbi Math", "Explanation of importance"),
                     tabPanel("Viterbi in Practice", "Explanation of importance"),
                     tabPanel("Trial Size Revisited", "Explanation of importance of simulation/sample size")
                     ))
      ## Other
      )
    )
  )

# Define server
server <- function(input, output) {
   
   output$dishonestCasinoBase <- renderPlot({
     set.seed(NULL)
     N = input$die_throws
     sim = HMM::simHMM(HMM_hmm, N) #simulate N dice throws using the HMM we made
     vit = HMM::viterbi(HMM_hmm, sim$observation) #take the observations, find the most likely sequence of hidden states that generated them
     f = HMM::forward(HMM_hmm, sim$observation)
     x = list(hmm = HMM_hmm, sim = sim, vit = vit)
     mn = "Fair vs unfair die (dishonest casino)"
     xlb = "Throw number"
     ylb = ""
     plot(x$sim$observation, ylim = c(-7.5, 6), pch = 3, main = mn,
          xlab = xlb, ylab = ylb, bty = "n", yaxt = "n")
     axis(2, at = 1:6)
     text(1, -1.4, adj = 0, cex = 0.8, col = "black", "True path: green = fair die, red = loaded die")
     for (i in 1:N) {
       if (x$sim$states[i] == "Fair")
         rect(i, -1, i + 1, 0, col = "#207d11", border = NA)
       else rect(i, -1, i + 1, 0, col = "#9e271e", border = NA)}
     text(1, -3.4, adj = 0, cex = 0.8, col = "black", "Most probable path")
     for (i in 1:N) {
       if (x$vit[i] == "Fair")
         rect(i, -3, i + 1, -2, col = "#207d11", border = NA)
       else rect(i, -3, i + 1, -2, col = "#9e271e", border = NA)
     }
     text(1, -5.4, adj = 0, cex = 0.8, col = "black", "Difference")
     differing = !(x$sim$states == x$vit)
     for (i in 1:N) {
       if (differing[i])
         rect(i, -5, i + 1, -4, col = rgb(0.3, 0.3, 0.3),
              border = NA)
       else rect(i, -5, i + 1, -4, col = rgb(0.9, 0.9, 0.9),
                 border = NA)
     }
     
     results = numeric(N)
     for (i in 1:N){
       ob<-x$sim$states[i]
       vit_pred<-x$vit[i]
       if (ob==vit_pred){
         results[i]<-1
       } else (results[i]<-0)
     }
     text(1, -7, adj = 0, cex = 0.8, col = "black", paste0("Total Agreement: ", round(100*mean(results),5),"%"), font=2)
   })
   
   
output$dishonestCasinoBoot <- renderPlot({
    N_throws<-input$die_throws_2
    N_trials<-input$casino_trials
    results<-numeric(N_trials)
    for (i in 1:N_trials){
      sim = HMM::simHMM(HMM_hmm, N_throws) #simulate N dice throws using the HMM we made
      vit = HMM::viterbi(HMM_hmm, sim$observation)
      df_x<-data.frame(1:N_throws)
      colnames(df_x)<-c("throw")
      df_x$sims<-sim$states
      df_x$viterbi<-vit
      df_x$value<-sim$observation
      differing = !(df_x$sims == df_x$viterbi)
      df_x$diff<-"No"
      df_x$diff[differing]<-"Yes"
      df_x<-df_x%>%
        mutate(num_diff=ifelse(grepl("No", diff), 0,1))
      results[i]<-df_x%>%
        dplyr::summarize(pct_same=round(mean(1-num_diff),3))%>%
        pull(pct_same)
    }
    data.frame(results)%>%
      mutate(trial=row_number())%>%
      ggplot()+
      geom_jitter(aes(x=trial, y=results, color=results), show.legend=F)+
      scale_y_continuous(labels = percent_format(), limits=c(0,1))+
      scale_color_steps()+
      labs(x="Trial", y="Accuracy")+
      annotate("label", x=N_trials/4, y=0.3, 
               label=paste0("Mean Accuracy: ", round(100*mean(results, na.rm=T), 3), "%"), fontface="bold",
               label.padding = unit(0.5, "lines")
      )+
      annotate("label", x=N_trials/4, y=0.2, 
               label=paste0("Range: ", round(100*min(results, na.rm=T), 3), "% - ", round(100*max(results, na.rm=T), 3), "%"), fontface="bold", label.padding = unit(0.5, "lines")
      )
  })


output$dishonestCasinoBootPlotly<-renderPlotly({
  N_throws<-input$die_throws_2
  N_trials<-input$casino_trials
  results<-numeric(N_trials)
  for (i in 1:N_trials){
    sim = HMM::simHMM(HMM_hmm, N_throws) #simulate N dice throws using the HMM we made
    vit = HMM::viterbi(HMM_hmm, sim$observation)
    df_x<-data.frame(1:N_throws)
    colnames(df_x)<-c("throw")
    df_x$sims<-sim$states
    df_x$viterbi<-vit
    df_x$value<-sim$observation
    differing = !(df_x$sims == df_x$viterbi)
    df_x$diff<-"No"
    df_x$diff[differing]<-"Yes"
    df_x<-df_x%>%
      mutate(num_diff=ifelse(grepl("No", diff), 0,1))
    results[i]<-df_x%>%
      dplyr::summarize(pct_same=round(mean(1-num_diff),3))%>%
      pull(pct_same)
  }
   
  test<-data.frame(results)%>%
    mutate(trial=row_number())
  mid<-0.5*N_trials
  plot_ly(test, x = ~trial, y = ~100*results, type = 'scatter',
          color = ~results, name="",
          hovertemplate = paste('<b>Trial Number: %{x} </b>',
                                '<br> Percent Correctly modeled: %{y}'))%>%
    layout(annotations = list(
      list(x = mid, y = 30, text = paste0("<b> Mean Accuracy: ", round(100*mean(results, na.rm=T), 2), "% </b>"), showarrow = F, 
           bgcolor="white", bordercolor="#0d0887", borderpad=10, borderwidth=3,
           font=list(size=12, color="black")), 
      list(x = mid, y = 20,
           text = paste0("<b> Range: ", round(100*min(results, na.rm=T), 2), "% - ", round(100*max(results, na.rm=T), 3), "% </b>"), showarrow = F, 
           font=list(size=10, color="black"))),
      yaxis = list(title = 'Percent Accurate', range=c(0, 105), ticksuffix="%"),
      xaxis = list(title="Trial Number")
    )
})

output$aphidBasePlot<-renderPlot({
  
  states <- c("Begin", "Fair", "Loaded")
  residues <- paste(1:6)
  ### Define transition probability matrix A
  A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9), nrow = 3)
  dimnames(A) <- list(from = states, to = states)
  ### Define emission probability matrix E
  E <- matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), nrow = 2, byrow = TRUE)
  dimnames(E) <- list(states = states[-1], residues = residues)
  ### Create the HMM object
  x <- structure(list(A = A, E = E), class = "HMM") #here we make an hmm, not as clean as HMM package
  ### Plot the model
  plot(x, textexp = 1.5)
  ### Optionally add the transition probabilities as text
  text(x = 0.02, y = 0.5, labels = "0.95")
  text(x = 0.51, y = 0.5, labels = "0.90")
  text(x = 0.5, y = 0.9, labels = "0.05")
  text(x = 0.5, y = 0.1, labels = "0.10")
  
})

output$woodmousePlot<-renderPlot({
  data(woodmouse)
  woodmouse.PHMM <- aphid::derivePHMM(woodmouse)
  plot(woodmouse.PHMM, from = 0, to = 5, main = "Example Profile HMM")
})
output$woodmouseData<-renderDataTable({
  data(woodmouse)
  x <- woodmouse[1:15, 1:10]
  as.character(cbind(x, fill.with.gaps = F))
  })
output$globinPlot<-renderPlot({
  data(globins)
  globins.PHMM <- aphid::derivePHMM(globins, residues = "AMINO", seqweights = NULL)
  plot(globins.PHMM, main = "Profile HMM for short globin alignment")
})
output$globinData<-renderDataTable({
  data(globins)
  print(globins)
})
output$globinsUnaligned<-renderText({
  data(globins)
  globins2 <- aphid::unalign(globins)
  print(globins2[1])
  #globins2[2]
})
output$globinsReAligned<-renderDataTable({
  data(globins)
  globins.PHMM <- aphid::derivePHMM(globins, residues = "AMINO", seqweights = NULL)
  globins2 <- aphid::unalign(globins)
  align(globins2, model = globins.PHMM, seqweights = NULL, residues = "AMINO")
})
output$firstThreeBirdsPHMM<-renderPlot({
  plot(birds.PHMM, from = 0, to = 10, main = "Three Birds profile HMM")
})
output$alignedBirds<-renderText({
  aligned<-align(birdsR, model=birds.PHMM)
  paste0("< Aligned an object of class ", class(aligned), " with sequences ",
         dimnames(aligned)[1], " and dimensions ", 
         dim(aligned)[1], "x", dim(aligned)[2], " >")
})
output$mvadHead<-renderDataTable({
  data("mvad", package = "TraMineR")
  head(mvad[,c(1,3,24:30)], n=3)
})
output$mvad_seq_plot<-renderPlot({
  TraMineR::seqplot(mvad_seq, with.legend=T)
  })
output$mavdHMMPlot<-renderPlot({
  plot(fit_hmm_mvad$model, combine.slices=0)
})
output$rdHMMPlot<-renderPlot({
  plot(fit_hmm_rd$model, combine.slices=0)
})
output$rdEmiss<-renderDataTable({
  rdemiss
})
output$rdTrans<-renderDataTable({
  rdtrans
})
output$depGGplot<-renderPlot({
  if(input$depSeedCheck==T){
  set.seed(input$depSeed)}
  #N_train<-as.numeric(input$depTraining)
  #N_test<-as.numeric(input$depTest)
  speed_train<-speed%>%
    slice_sample(n=input$depTraining, replace=T)
  speed_test<-speed%>%
    slice_sample(n=input$depTest, replace=T)
  xDepth<-0.75*input$depTest
  mod <- depmix(response = rt ~ 1, #response to be modelled.
                data = speed_train, #optional data to use, a data frame 
                nstates = 2)#number of hidden states

  fmspeed <- fit(mod, emcontrol=em.control(rand=FALSE))
  
  modNew <- depmix(rt~1, data=speed_test, nstates=2) #now make a new model of the test data
  modNew <- setpars(modNew,getpars(fmspeed)) #use the parameters from the trained model
  modNew <- fit(modNew)#fit it
  
  predStates <- depmixS4::posterior(modNew) # predictions of the hidden states
  
  speed_test<-speed_test%>%
    mutate(true_state=ifelse(grepl("inc", corr), 2, 1))
  check<-cbind(speed_test, predStates)
  check<-check%>%
    dplyr::rename(pred_state=state)%>%
    mutate(correct=ifelse(true_state==pred_state, "True", "False"))%>%
    mutate(Prediction_number=row_number())%>%
    mutate(agreed=ifelse(correct=="True", 1, 0))%>%
    dplyr::select(rt, corr, true_state, pred_state, correct, Prediction_number, agreed)
  
  check%>%
    ggplot(aes(x=Prediction_number, y=correct, color=correct))+
    geom_point(shape=3, size=3, show.legend=F)+
    scale_color_manual(values=c("#9e271e", "#207d11"))+
    labs(x="Prediction Number", y="")+
    ggtitle("Speed Data HMM Results")+
    annotate("label", x=xDepth, y=1.75, 
             label=paste0("Correct: ", round(100*mean(check$agreed, na.rm=T), 2), "%"), fontface="bold", label.padding = unit(0.5, "lines")
    )+
    annotate("label", x=xDepth, y=0.75, 
             label=paste0("Incorrect: ", round(100*(1-mean(check$agreed, na.rm=T)), 2), "%"), fontface="bold", label.padding = unit(0.5, "lines")
    )
  
})
output$depPlotly<-renderPlot({
  if(input$depSeedCheck==T){
    set.seed(input$depSeed)}
  N_train<-as.numeric(input$depTraining)
  N_test<-as.numeric(input$depTest)
  speed_train<-speed%>%
    slice_sample(n=N_train, replace=T)
  speed_test<-speed%>%
    slice_sample(n=N_test, replace=T)
  xDepth<-0.75*N_test
  mod <- depmix(response = rt ~ 1, #response to be modelled.
                data = speed_train, #optional data to use, a data frame 
                nstates = 2)#number of hidden states
  
  fmspeed <- fit(mod, emcontrol=em.control(rand=FALSE))
  
  modNew <- depmix(rt~1, data=speed_test, nstates=2) #now make a new model of the test data
  modNew <- setpars(modNew,getpars(fmspeed)) #use the parameters from the trained model
  modNew <- fit(modNew)#fit it
  
  predStates <- depmixS4::posterior(modNew) # predictions of the hidden states
  
  speed_test<-speed_test%>%
    mutate(true_state=ifelse(grepl("inc", corr), 2, 1))
  check<-cbind(st2, predStates)
  check<-check%>%
    dplyr::rename(pred_state=state)%>%
    mutate(correct=ifelse(true_state==pred_state, "True", "False"))%>%
    mutate(Prediction_number=row_number())%>%
    mutate(agreed=ifelse(correct=="True", 1, 0))%>%
    dplyr::select(rt, corr, true_state, pred_state, correct, Prediction_number, agreed)
  
  
  fig<-check%>%
    ggplot(aes(x=Prediction_number, y=correct, color=correct))+
    geom_point(shape=3, size=3, show.legend=F)+
    scale_color_manual(values=c("#9e271e", "#207d11"))+
    labs(x="Prediction Number", y="")+
    ggtitle("Speed Data HMM Results")
  fig<-ggplotly(fig)
  
  fig%>%
    layout(annotations = list(
      list(x = xDepth, y = 1.75, text = paste0("<b> Percentage Correct: ", round(100*mean(check$agreed, na.rm=T), 2), "% </b>"), showarrow = F, 
           bgcolor="white", bordercolor="#207d11", borderpad=10, borderwidth=3,
           font=list(size=12, color="black")),
      list(x = xDepth, y = 0.75, text = paste0("<b> Percentage Incorrect: ", round(100*(1-mean(check$agreed, na.rm=T)), 2), "% </b>"), showarrow = F, 
           bgcolor="white", bordercolor="#9e271e", borderpad=10, borderwidth=3,
           font=list(size=12, color="black"))
    ),
    yaxis = list(title = ' ', range=c(0, 2.5), ticksuffix="%")
    )
})
output$sp500DT<-renderDataTable({
  head(sp500, 3)
})
output$sp500TSClose<-renderPlot({
  plot(ts(sp500$Close, start=c(1950,1), deltat=1/12), main="Closing Prices", ylab="Price")
})
output$sp500TSLR<-renderPlot({
  plot(ts(sp500$logret, start=c(1950,1), deltat=1/12), main="Log Return", ylab="Log Return")
})
output$poissonBIC<-renderPlot({
  set.seed(input$bicSeedCheck)
  y1 <- rpois(n=input$length1, lambda=input$lambda1)
  y2 <- rpois(n=input$length2, lambda=input$lambda2)
  ydf <- data.frame(y=c(y1,y2))
  
  m1 <- depmix(y~1, nstates=1,family=poisson(), data=ydf)
  set.seed(1)
  fm1 <- fit(m1)
  m2 <- depmix(y~1, nstates=2,family=poisson(), data=ydf)
  set.seed(1)
  fm2 <- fit(m2)
  m3 <- depmix(y~1, nstates=3,family=poisson(), data=ydf)
  set.seed(1)
  fm3 <- fit(m3,em=em.control(maxit=200))
  m4 <- depmix(y~1, nstates=4,family=poisson(), data=ydf)
  set.seed(1)
  fm4 <- fit(m4,em=em.control(maxit=200))
  m5 <- depmix(y~1, nstates=5,family=poisson(), data=ydf)
  set.seed(1)
  fm5 <- fit(m5,em=em.control(maxit=200))
  
  # plot the BICs to select the proper model
  plot(1:5, c(BIC(fm1),BIC(fm2),BIC(fm3),BIC(fm4),BIC(fm5)),
       ty="b", ylab="BIC", xlab="Number of Hidden States")
})

} #close server
# Run the application 
shinyApp(ui = ui, server = server)

