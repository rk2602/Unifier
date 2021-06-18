#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Contact Nicholas Parham (NP) at nick-99@att.net for comments or corrections.

library(shiny)
library(readxl)
library(dplyr)
library(DT)
library(sjmisc)
library(PeriodicTable)
library(data.table)

   #########################
###### FILE DEPENDENCIES ######
   #########################

elements.list = read.csv('elementlist.csv') # should be a dependency
elements = elements.list$symbols

   #################
###### FUNCTIONS ######
   #################

s_element = function(element){ # isolate element name for use in mass() function
    if (unlist(stringr::str_split(element, '[(]'), use.names = F)[1] %in% elements){
        if (length(unlist(stringr::str_split(element, '[+]'), use.names = F)) > 1){
            symbol = unlist(stringr::str_split(element, '[(]'), use.names = F)[1]
        }else if (length(unlist(stringr::str_split(element, '[-]'), use.names = F)) > 1){
            symbol = unlist(stringr::str_split(element, '[(]'), use.names = F)[1]
        }else{
            symbol = 'H' # default if element not found
        }
        if (is.null(symbol) | is.na(symbol)){
            symbol = 'H'
        }
    }else{
        symbol = 'H'
    }
    return(symbol)
}

   ####################
###### SERVER LOGIC ######
   ####################

shinyServer(function(input, output) {
    
       ################
    ###### HOME TAB ######
       ################
    
    missingInput = eventReactive(input$scan, {
        filepath = input$dataset.test
        dataset = read.csv(filepath$datapath)
        filepath2 = input$mineral.test
        mineral.ref = read_excel(filepath2$datapath)
        
        minerals = unique(c(dataset$Mineral, 
                            dataset$Electrolyte1,
                            dataset$Electrolyte2,
                            dataset$Electrolyte3,
                            dataset$Electrolyte4,
                            dataset$Electrolyte5,
                            dataset$Electrolyte6,
                            dataset$Sorbent))
        ref.minerals = unique(mineral.ref$minerals)
        elements.list = read.csv('elementlist.csv') # should be a dependency
        elements = elements.list$symbols
        
        # filter out elements
        for (i in c(1:length(minerals))){
            if (unlist(stringr::str_split(minerals[i], '[(]'), use.names = F)[1] %in% elements){
                if (length(unlist(stringr::str_split(minerals[i], '[+]'), use.names = F)) > 1){
                    minerals[i] = 'Element'
                }else if (length(unlist(stringr::str_split(minerals[i], '[-]'), use.names = F)) > 1){
                    minerals[i] = 'Element'
                }else{
                    # not element
                }
            }
        }

        missing = minerals[!(minerals %in% ref.minerals) & !(minerals == '') & !(is.na(minerals)) 
                           & !(minerals == 'Element') & !(minerals == 'pH')]
        missing = data.frame(missing = missing) # output missing from mineral-ref
        missing
    })
    
    output$missing = renderTable({
        missing = missingInput()
        missing
    })

       ####################
    ###### UNIFIER TAB #######
       ####################
    
    sc.datasetInput = reactive({ # read in Dataset.xlsx as dataframe
        filepath = input$sc.dataset
        sc.dataset = read.csv(filepath$datapath)
        sc.dataset
    })
    
    sc.dataInput = reactive({ # read in Data.xlsx as dataframe
        filepath = input$sc.data
        sc.data = read.csv(filepath$datapath)
        sc.data
    })
    
    sc.mineralsInput = reactive({ # read in mineral-ref.xlsx as dataframe
        filepath = input$sc.minerals
        sc.minerals = read_excel(filepath$datapath, 1)
        sc.minerals
    })
    
    sc.dataset = eventReactive(input$unify, { # create uniform dataset
        
        #########################
        ### READ IN USER SD's ###
        #########################
        sd1 = input$sd1       # temp SD
        sd2 = input$sd2 / 100 # elyte linear SD
        sd3 = input$sd3 / 100 # elyte log SD
        sd4 = input$sd4       # pH SD
        sd5 = input$sd5 / 100 # mineral SD
        sd6 = input$sd6 / 100 # mineralSA SD
        sd7 = input$sd7 / 100# mineral sites SD
        sd8 = input$sd8 / 100 # CEC SD
        sd9 = input$sd9 / 100 # gas SD
        sd10 = input$sd10 / 100 # sorbent SD
        sd11 = input$sd11 / 100 # charge linear SD
        sd12 = input$sd12 / 100# charge log SD
        sd13 = input$sd13 / 100 # sorbed %/frac SD
        sd14 = input$sd14 / 100 # sorbed Kd/Rd linear SD
        sd15 = input$sd15 / 100 # sorbed Kd/Rd log SD
        sd16 = input$sd16 / 100 # sorbed linear SD
        sd17 = input$sd17 / 100 # sorbed log SD
        
        ###############################
        ### READ IN USER DATA FILES ###
        ###############################
        dataset = sc.datasetInput()
        dat = sc.dataInput()
        mineral.ref = sc.mineralsInput()
        dataset = right_join(dataset, dat, by = 'Set') # join Dataset and Data on Set column
        dataset = dataset[dataset$Mineral != '' & !is.na(dataset$Mineral),] # remove incomplete data
        p = nrow(dataset) # get rows to create empty dataframes and control loops
        
        for (i in c(1:ncol(dataset))){ # fix formatting assumption of logical for numeric columns
            if (is.logical(dataset[,i])){
                dataset[,i] = as.numeric(dataset[,i])
            }
        }
        
        ########################################
        ### OUTPUT DATAFRAME STRUCTURE SETUP ###
        ########################################
        refs = dataset$Reference # transfer without conversion
        sets = dataset$Set # transfer without conversion
        setIDs = dataset$number # transfer without conversion
        minerals = dataset$Mineral # transfer without conversion
        formulas = dataset$Mineral_Formula # transfer without conversion
        sources = dataset$Mineral_source # transfer without conversion
        temps = dataset$Temp # transfer without conversion
        sd.temps = c(rep(NA,p))
        
        ###############################
        ### HANDLE TEMP ESTIMATIONS ###
        ###############################
        
        for (i in c(1:p)){
            if (temps[i] != '' & !is.na(temps[i])){
                if (!is.na(sd.temps[i])){
                    # good
                }else{
                    new.sd = sd1 # 7JUL20
                    sd.temps[i] = new.sd
                }
            }
        }
        
        ##################################################################
        ### HANDLE ELECTROLYTE CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        ##################################################################
        electrolytes = data.frame(Electrolyte1 = rep(NA,p), Electrolyte1_val = rep(NA,p), 
                                  Electrolyte1_SD = rep(NA,p), Electrolyte1_units = rep(NA,p),
                                  Electrolyte2 = rep(NA,p), Electrolyte2_val = rep(NA,p),
                                  Electrolyte2_SD = rep(NA,p), Electrolyte2_units = rep(NA,p),
                                  Electrolyte3 = rep(NA,p), Electrolyte3_val = rep(NA,p),
                                  Electrolyte3_SD = rep(NA,p), Electrolyte3_units = rep(NA,p),
                                  Electrolyte4 = rep(NA,p), Electrolyte4_val = rep(NA,p),
                                  Electrolyte4_SD = rep(NA,p), Electrolyte4_units = rep(NA,p),
                                  Electrolyte5 = rep(NA,p), Electrolyte5_val = rep(NA,p),
                                  Electrolyte5_SD = rep(NA,p), Electrolyte5_units = rep(NA,p),
                                  Electrolyte6 = rep(NA,p), Electrolyte6_val = rep(NA,p),
                                  Electrolyte6_SD = rep(NA,p), Electrolyte6_units = rep(NA,p),
                                  Electrolyte7 = rep(NA,p), Electrolyte7_val = rep(NA,p),
                                  Electrolyte7_SD = rep(NA,p), Electrolyte7_units = rep(NA,p),
                                  pH = rep(NA,p), pH_SD = rep(NA,p)) # create electrolytes sub-dataframe for future use
        
        pHs = electrolytes$pH # read values
        sd.pHs = electrolytes$pH_SD # read values
        
        for (n in c(1:7)){ # perform electrolyte conversions, pH transfers, and SD estimations
            elyte = paste('Electrolyte', n, sep = '')
            elyte.val = paste('Electrolyte', n, '_val', sep = '')
            elyte.sd = paste('Electrolyte', n, '_SD', sep = '')
            elyte.unit = paste('Electrolyte', n, '_units', sep = '')
            elytes = as.vector(dataset[,elyte]) # read values
            val.elytes = as.vector(dataset[,elyte.val]) # read values
            sd.elytes = as.vector(dataset[,elyte.sd]) # read values
            units.elytes = as.vector(dataset[,elyte.unit]) # read values
            
            for (i in c(1:p)){
                if (units.elytes[i] != '' & !is.na(units.elytes[i])){
                    if (units.elytes[i] == 'pH'){
                        pHs[i] = val.elytes[i]
                        if (!is.na(sd.elytes[i])){
                            sd.pHs[i] = sd.elytes[i]
                        }else{
                            sd.pHs[i] = sd4 # 7JUL20
                        }
                        elytes[i] = NA
                        val.elytes[i] = NA
                        sd.elytes[i] = NA
                        units.elytes[i] = NA
                        
                    }else if (units.elytes[i] == 'mg/L'){
                        if (elytes[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == elytes[i], 'masses'])
                            new.val = val.elytes[i] / (1000 * molar.mass)
                            if (!is.na(sd.elytes[i])){
                                new.sd = sd.elytes[i] / (1000 * molar.mass)
                                sd.elytes[i] = new.sd
                            }else{
                                new.sd = new.val * sd2 # 7JUL20
                                sd.elytes[i] = new.sd
                            }
                            val.elytes[i] = new.val
                            units.elytes[i] = 'mol/L'
                            
                        }else{
                            if (elytes[i] == 'Zn(+2)'){
                                molar.mass = 65.38
                            }else{
                                molar.mass = mass(s_element(elytes[i]))
                            }
                            new.val = (val.elytes[i]) / (1000 * molar.mass)
                            if (!is.na(sd.elytes[i])){
                                new.sd = (sd.elytes[i]) / (1000 * molar.mass)
                                sd.elytes[i] = new.sd
                            }else{
                                new.sd = new.val * sd2 # 7JUL20
                                sd.elytes[i] = new.sd
                            }
                            val.elytes[i] = new.val
                            units.elytes[i] = 'mol/L'
                        }
                        
                    }else if (units.elytes[i] == 'mmol/L'){
                        new.val = val.elytes[i] / 1000
                        if (!is.na(sd.elytes[i])){
                            new.sd = sd.elytes[i] / 1000
                            sd.elytes[i] = new.sd
                        }else{
                            new.sd = new.val * sd2 # 7JUL20
                            sd.elytes[i] = new.sd
                        }
                        val.elytes[i] = new.val
                        units.elytes[i] = 'mol/L'
                        
                    }else if (units.elytes[i] == 'mol/l'){
                        if (!is.na(sd.elytes[i])){
                            # good
                        }else{
                            new.sd = val.elytes[i] * sd2 # 7JUL20
                            sd.elytes[i] = new.sd
                        }
                        units.elytes[i] = 'mol/L'
                        
                    }else if (units.elytes[i] == 'mol/kg'){
                        if (!is.na(sd.elytes[i])){
                            # good
                        }else{
                            new.sd = val.elytes[i] * sd2 # 7JUL20
                            sd.elytes[i] = new.sd
                        }
                        units.elytes[i] = 'mol/L'
                        
                    }else if (units.elytes[i] == 'mol/L'){
                        # mol/L good
                        if (!is.na(sd.elytes[i])){
                            # good
                        }else{
                            new.sd = val.elytes[i] * sd2 # 7JUL20
                            sd.elytes[i] = new.sd
                        }
                    }
                }
            }
            electrolytes[,elyte] = elytes # update values
            electrolytes[,elyte.val] = val.elytes # update values
            electrolytes[,elyte.sd] = sd.elytes # update values
            electrolytes[,elyte.unit] = units.elytes # update values
        }
        
        electrolytes$pH = pHs # update values
        electrolytes$pH_SD = sd.pHs # update values
        
        ##############################################################
        ### HANDLE SORBENT CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        ##############################################################
        sorbents = dataset$Sorbent  # read values
        val.sorbents = dataset$Sorbent_val  # read values
        sd.sorbents = dataset$Sorbent_SD  # read values
        units.sorbents = dataset$sorbent.units  # read values
        
        for (i in c(1:p)){ # perform sorbent conversions and SD estimations
            if (units.sorbents[i] != '' & !is.na(units.sorbents[i])){
                if (sorbents[i] %in% mineral.ref$minerals){
                    molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == sorbents[i], 'masses'])
                }else{
                    molar.mass = mass(s_element(sorbents[i]))
                }
                if (units.sorbents[i] %in% c('ppb', 'ug/L')){
                    new.val = (val.sorbents[i] * 0.001) / (1000 * molar.mass)
                    if (!is.na(sd.sorbents[i])){
                        new.sd = (sd.sorbents[i] * 0.001) / (1000 * molar.mass)
                        sd.sorbents[i] = new.sd
                    }else{
                        new.sd = new.val * sd10 # 7JUL20
                        sd.sorbents[i] = new.sd
                    }
                    val.sorbents[i] = new.val
                    units.sorbents[i] = 'mol/L'
                    
                }else if (units.sorbents[i] %in% c('ppm', 'mg/L')){
                    new.val = (val.sorbents[i]) / (1000 * molar.mass)
                    if (!is.na(sd.sorbents[i])){
                        new.sd = (sd.sorbents[i]) / (1000 * molar.mass)
                        sd.sorbents[i] = new.sd
                    }else{
                        new.sd = new.val * sd10 # 7JUL20
                        sd.sorbents[i] = new.sd
                    }
                    val.sorbents[i] = new.val
                    units.sorbents[i] = 'mol/L'
                    
                }else if (units.sorbents[i] == 'mol/kg'){
                    if (!is.na(sd.sorbents[i])){
                        # good
                    }else{
                        new.sd = val.sorbents[i] * sd10 # 7JUL20
                        sd.sorbents[i] = new.sd
                    }
                    units.sorbents[i] = 'mol/L'
                    
                }else if (units.sorbents[i] == 'mol/L'){
                    # mol/L good
                    if (!is.na(sd.sorbents[i])){
                        # good
                    }else{
                        new.sd = val.sorbents[i] * sd10 # 7JUL20
                        sd.sorbents[i] = new.sd
                    }
                    
                }else if (units.sorbents[i] == 'mmol/L'){
                    new.val = val.sorbents[i] / 1000
                    if (!is.na(sd.sorbents[i])){
                        new.sd = sd.sorbents[i] / 1000
                        sd.sorbents[i] = new.sd
                    }else{
                        new.sd = new.val * sd10 # 7JUL20
                        sd.sorbents[i] = new.sd
                    }
                    val.sorbents[i] = new.val
                    units.sorbents[i] = 'mol/L'
                }
            }
        }
        
        sorbent = data.frame(Sorbent = sorbents, Sorbent_val = val.sorbents, Sorbent_SD = sd.sorbents, 
                             Sorbent_units = units.sorbents) # update values and create sub-dataframe for future use
        
        aqueous = data.frame(Aq_val = rep(NA,p), Aq_SD = rep(NA,p), 
                             Aq_units = rep(NA,p)) # create aqueous sub-dataframe for future use
        sorbed = data.frame(Sorbed_val = rep(NA,p), Sorbed_SD = rep(NA,p), 
                            Sorbed_units = rep(NA,p)) # create sorbed sub-dataframe for future use
        
        ##############################################################
        ### HANDLE MINERAL CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        ##############################################################
        val.minerals = dataset$Mineral_val # read values
        sd.minerals = dataset$Mineral_SD # read values
        units.minerals = dataset$Mineral_units # read values
        
        for (i in c(1:p)){
            if (units.minerals[i] != '' & !is.na(units.minerals[i])){ # perform mineral conversions and SD estimations
                if (units.minerals[i] %in% c('g/l','g/L')){
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        new.val = val.minerals[i] / molar.mass
                        if (!is.na(sd.minerals[i])){
                            new.sd = sd.minerals[i] / molar.mass
                            sd.minerals[i] = new.sd
                        }else{
                            new.sd = new.val * sd5 # 7JUL20
                            sd.minerals[i] = new.sd
                        }
                        val.minerals[i] = new.val
                        units.minerals[i] = 'mol/L'
                    }
                }
            }
        }
        
        mineral = data.frame(Mineral_val = val.minerals, Mineral_SD = sd.minerals, 
                             Mineral_units = units.minerals) # update values and create sub-dataframe for future use
        
        #################################################################
        ### HANDLE MINERAL SA CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        #################################################################
        val.mineralSA = mineralSA = dataset$MineralSA # read values
        sd.mineralSA = dataset$MineralSA_SD # read values
        units.mineralSA = dataset$MineralSA_units # read values
        
        for (i in c(1:p)){
            if (units.mineralSA[i] != '' & !is.na(units.mineralSA[i])){ # perform mineral SA conversions and SD estimations
                if (units.mineralSA[i] == 'mm'){
                    if(minerals[i] %in% mineral.ref$minerals){
                        d = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'densities'])
                        new.val = ((6 * val.mineralSA[i]^2) / 1000000 ) / ((val.mineralSA[i]^3 * d) / (1000))
                        if (!is.na(sd.mineralSA[i])){
                            new.sd = ((6 * sd.mineralSA[i]^2) / 1000000 ) / ((sd.mineralSA[i]^3 * d) / (1000))
                            sd.mineralSA[i] = new.sd
                        }else{
                            new.sd = new.val * sd6 # 7JUL20
                            sd.mineralSA[i] = new.sd
                        }
                        val.mineralSA[i] = new.val
                        units.mineralSA[i] = 'm2/g'
                    }
                    
                }else if (units.mineralSA[i] == 'm2/g'){
                    # m2/g good
                    if (!is.na(sd.mineralSA[i])){
                        # good
                    }else{
                        new.sd = val.mineralSA[i] * sd6 # 7JUL20
                        sd.mineralSA[i] = new.sd
                    }
                }
            }
        }
        
        mineralSA = data.frame(MineralSA = val.mineralSA, MineralSA_SD = sd.mineralSA, 
                               MineralSA_units = units.mineralSA) # update values and create sub-dataframe for future use
        
        sites = dataset[,c('Mineralsites', 'Mineralsites_SD', 'Mineralsites_units')] # transfer without conversion
        cec = dataset[,c('CEC', 'CEC_SD', 'CEC_units')] # transfer without conversion
        
        gases = dataset[,c('Gas1', 'Gas1_val', 'Gas1_SD', 'Gas1_units', # transfer without conversion
                           'Gas2', 'Gas2_val', 'Gas2_SD', 'Gas2_units', # transfer without conversion
                           'Gas3', 'Gas3_val', 'Gas3_SD', 'Gas3_units')] # transfer without conversion
        
        charges = data.frame(SurfCharge_val = rep(NA,p), SurfCharge_SD = rep(NA,p), SurfCharge_units = rep(NA,p)) 
        
        ###########################
        ### HANDLE SITES VALUES ###
        ###########################
        val.sites = dataset$Mineralsites # read values
        
        # insert code based off of user input
        if (input$sites == 'Replace'){
            for (i in c(1:p)){
                if (minerals[i] %in% mineral.ref$minerals){
                    val.sites[i] = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'sites'])
                }
            }
        }else if (input$sites == 'Fill'){
            for (i in c(1:p)){
                if (val.sites[i] != '' & !is.na(val.sites[i])){
                    # good
                }else{
                    if (minerals[i] %in% mineral.ref$minerals){
                        val.sites[i] = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'sites'])
                    }
                }
            }
        }
        
        dataset$Mineralsites = val.sites # update values
        sites$Mineralsites = val.sites # update values
        
        ##############################
        ### HANDLE CEC ESTIMATIONS ###
        ##############################
        val.cec = dataset$CEC # read values
        sd.cec = dataset$CEC_SD # read values
        
        for (i in c(1:p)){
            if (val.cec[i] != '' & !is.na(val.cec[i])){
                if (!is.na(sd.cec[i])){
                    # good
                }else{
                    new.sd = val.cec[i] * sd8 # 7JUL20
                    sd.cec[i] = new.sd
                }
            }
        }
        
        cec$CEC_SD = sd.cec # update values
        
        ########################################
        ### HANDLE MINERAL SITES ESTIMATIONS ###
        ########################################
        val.sites = dataset$Mineralsites # read values
        sd.sites = dataset$Mineralsites_SD # read values
        units.sites = dataset$Mineralsites_units # read values
        
        for (i in c(1:p)){
            if (val.sites[i] != '' & !is.na(val.sites[i])){
                if (!is.na(sd.sites[i])){
                    # good
                }else{
                    new.sd = val.sites[i] * sd7 # 7JUL20
                    sd.sites[i] = new.sd
                }
                units.sites[i] = 'sites/nm2'
            }
        }
        
        sites$Mineralsites_SD = sd.sites # update values
        sites$Mineralsites_units = units.sites # update values
        
        ##############################
        ### HANDLE GAS ESTIMATIONS ###
        ##############################
        for (n in c(1:3)){
            gas.val = paste('Gas', n, '_val', sep = '')
            gas.sd = paste('Gas', n, '_SD', sep = '')
            val.gases = as.vector(gases[,gas.val]) # read values
            sd.gases = as.vector(gases[,gas.sd]) # read values
            
            for (i in c(1:p)){
                if (!is.na(sd.gases[i])){
                    # good
                }else{
                    new.sd = val.gases[i] * sd9 # 7JUL20
                    sd.gases[i] = new.sd
                }
            }
            gases[,gas.sd] = sd.gases # update values
        }
        
        ############################################
        ### BIND AND MERGE INTO OUTPUT STRUCTURE ###
        ############################################
        a = data.frame(Reference = refs, Set = sets, SetID = setIDs, Mineral = minerals, Mineral_formula = formulas, 
                       Mineral_source = sources, Temp = temps, Temp_SD = sd.temps)
        b = cbind(a, electrolytes, sorbent, aqueous, sorbed, charges, mineral, mineralSA, sites, cec, gases)
        c = cbind(b, dataset[,c('X_axis', 'X_val', 'X_SD', 'X_units', 'Y_axis', 'Y_val', 'Y_SD', 'Y_units')])
        dataset = c # re-define dataset as output
        
        ##################################################
        ### HANDLE X/Y-AXIS pH CONVERSIONS & TRANSFERS ###
        ##################################################
        axis.y = dataset$Y_axis # read values
        val.y = dataset$Y_val # read values
        sd.y = dataset$Y_SD # read values
        units.y = dataset$Y_units # read values
        
        axis.x = dataset$X_axis # read values
        val.x = dataset$X_val # read values
        sd.x = dataset$X_SD # read values
        units.x = dataset$X_units # read values
        
        pHs = dataset$pH # read values
        sd.pHs = dataset$pH_SD # read values
        
        for (i in c(1:p)){
            if (!is.na(axis.y[i])){
                if (axis.y[i] == 'pH'){
                    pHs[i] = val.y[i]
                    if (!is.na(sd.y[i])){
                        sd.pHs[i] = sd.y[i]
                    }else{
                        sd.pHs[i] = sd4 # 7JUL20
                    }
                    axis.y[i] = NA
                    val.y[i] = NA
                    sd.y[i] = NA
                    units.y[i] = NA
                }
            }
            
            if (!is.na(axis.x[i])){
                if (axis.x[i] == 'pH' | axis.x[i] == '-log(mol/L)' | units.x[i] == '-log(mol/L)'){ # 19JUL20
                    pHs[i] = val.x[i]
                    if (!is.na(sd.x[i])){
                        sd.pHs[i] = sd.x[i]
                    }else{
                        sd.pHs[i] = sd4 # 7JUL20
                    }
                    axis.x[i] = NA
                    val.x[i] = NA
                    sd.x[i] = NA
                    units.x[i] = NA
                }
            }
        }
        
        dataset$Y_axis = axis.y # update values
        dataset$Y_val = val.y # update values
        dataset$Y_SD = sd.y  # update values
        dataset$Y_units = units.y # update values
        
        dataset$X_axis = axis.x # update values
        dataset$X_val = val.x # update values
        dataset$X_SD = sd.x # update values
        dataset$X_units = units.x # update values
        
        dataset$pH = pHs
        dataset$pH_SD = sd.pHs
        
        #############################################################
        ### HANDLE Y-AXIS CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        #############################################################
        axis.y = dataset$Y_axis # read values
        val.y = dataset$Y_val # read values
        sd.y = dataset$Y_SD # read values
        units.y = dataset$Y_units # read values
        
        axis.x = dataset$X_axis # read values
        val.x = dataset$X_val # read values
        sd.x = dataset$X_SD # read values
        units.x = dataset$X_units # read values
        
        areas = dataset$MineralSA # read values
        sorbents = dataset$Sorbent # read values
        val.minerals = dataset$Mineral_val # read values
        
        pHs = dataset$pH # read values
        sd.pHs = dataset$pH_SD # read values
        
        val.charges = dataset$SurfCharge_val # read values
        sd.charges = dataset$SurfCharge_SD # read values
        units.charges = dataset$SurfCharge_units # read values
        
        # clean y data values
        for (i in c(1:p)){
            if (units.y[i] != '' & !is.na(units.y[i])){
                if (units.y[i] == 'log(L/kg)'){
                    new.val = 10^(val.y[i])
                    if (!is.na(sd.y[i])){
                        new.sd = 10^(sd.y[i])
                        sd.y[i] = new.sd
                    }
                    val.y[i] = new.val
                    axis.y[i] = 'log(Rd)'
                    units.y[i] = 'L/kg'
                    
                }else if (units.y[i] == 'log(mL/g)' | units.y[i] == 'log(ml/g)'){
                    new.val = 10^(val.y[i])
                    if (!is.na(sd.y[i])){
                        new.sd = 10^(sd.y[i])
                        sd.y[i] = new.sd
                    }
                    val.y[i] = new.val
                    axis.y[i] = 'log(Kd)'
                    units.y[i] = 'mL/g'
                    
                }else if (units.y[i] == 'log(mol/L)'){
                    new.val = 10^(val.y[i])
                    if (!is.na(sd.y[i])){
                        new.sd = 10^(sd.y[i])
                        sd.y[i] = new.sd
                    }else{
                        new.sd = new.val * sd17 # 7JUL20
                        sd.y[i] = new.sd
                    }
                    val.y[i] = new.val
                    axis.y[i] = 'aqueous'
                    units.y[i] = 'mol/L'
                    
                }else if (units.y[i] == 'log(mol/kg)' | units.y[i] == 'log(mol/m2)' | units.y[i] == 'log(mol/g)'){
                    new.val = 10^(val.y[i])
                    if (!is.na(sd.y[i])){
                        new.sd = 10^(sd.y[i])
                        sd.y[i] = new.sd
                    }else{
                        new.sd = new.val * sd17 # 7JUL20
                        sd.y[i] = new.sd
                    }
                    val.y[i] = new.val
                    axis.y[i] = 'sorbed'
                    if (units.y[i] == 'log(mol/kg)'){
                        units.y[i] = 'mol/kg'
                        
                    }else if (units.y[i] == 'log(mol/m2)'){
                        units.y[i] = 'mol/m2'
                        
                    }else if (units.y[i] == 'log(mol/g)'){
                        units.y[i] = 'mol/g'
                        
                    }
                    
                }else if (units.y[i] == 'log(mol/mol)'){
                    new.val = 10^(val.y[i])
                    if (!is.na(sd.y[i])){
                        new.sd = 10^(sd.y[i])
                        sd.y[i] = new.sd
                    }else{
                        new.sd = new.val * sd17 # 7JUL20
                        sd.y[i] = new.sd
                    }
                    val.y[i] = new.val
                    axis.y[i] = 'sorbed'
                    units.y[i] = 'mol/mol'
                    
                    ### SURFACE CHARGE CONVERSIONS AND TRANSFERS ###
                }else if (units.y[i] == 'C/m2'){
                    new.val = val.y[i] * 1000000 * (1 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * 1000000 * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'microC/m2'){
                    new.val = val.y[i] * (1 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'C/g'){
                    new.val = val.y[i] * (1 / areas[i]) * 1000000 * (1 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * (1 / areas[i]) * 1000000 * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'eq/g'){
                    new.val = val.y[i] * 96485 * (1 / areas[i]) * 1000000 * (1 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * 96485 * (1 / areas[i]) * 1000000 * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'meq/kg'){
                    new.val = val.y[i] * 96485 * (1 / areas[i]) * (1 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * 96485 * (1 / areas[i]) * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'cmol/kg'){
                    new.val = val.y[i] * (1 / 100) * (1 / 1000) * 96485 * (1 / areas[i]) * (1000000 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * (1 / 100) * (1 / 1000) * 96485 * (1 / areas[i]) * (1000000 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'micromol/m2' & axis.y[i] != 'sorbed'){ # 19JUL20
                    new.val = val.y[i] * 96485 * (1 / 10000) # 19JUL20
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * 96485 * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = new.val * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (axis.y[i] == 'total H(+)'){
                    if (units.y[i] == 'mmol/L'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        new.val = ((val.y[i]/1000) - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                        if (!is.na(sd.y[i])){
                            new.sd = ((sd.y[i]/1000) - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 7JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                        
                    }else if (units.y[i] == 'mol/L'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        new.val = (val.y[i] - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                        if (!is.na(sd.y[i])){
                            new.sd = (sd.y[i] - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 7JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                    }
                    
                }else if (units.y[i] == 'microC/m2'){
                    new.val = val.y[i] * (1 / 10000)
                    if (!is.na(sd.y[i])){
                        new.sd = sd.y[i] * (1 / 10000)
                        sd.charges[i] = new.sd
                    }else{
                        new.sd = val.y[i] * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    
                }else if (units.y[i] == 'microC/cm2'){
                    if (!is.na(sd.y[i])){
                        # good
                    }else{
                        new.sd = val.y[i] * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = val.y[i]
                    units.charges[i] = 'microC/cm2'
                    
                }else if (axis.y[i] == 'charge'){
                    if (units.y[i] == 'mol/L'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        new.val = val.y[i] * 96485 * (1 / min.val) * (1 / areas[i]) * (1000000 / 10000)
                        if (!is.na(sd.y[i])){
                            new.sd = sd.y[i] * 96485 * (1 / min.val) * (1 / areas[i]) * (1000000 / 10000)
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 16JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                        
                    }else if (units.y[i] == 'mol/kg'){ # 16JUL20
                        new.val = val.y[i] * 96485 * (1 / 1000) * (1 / areas[i]) * (1000000 / 10000)
                        if (!is.na(sd.y[i])){
                            new.sd = sd.y[i] * 96485 * (1 / 1000) * (1 / areas[i]) * (1000000 / 10000)
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 7JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                    }
                }
                
                if (axis.y[i] == 'adsorbed'){
                    axis.y[i] = 'sorbed'
                }
                
                # sorbed prep
                if (axis.y[i] == 'sorbed'){
                    if (units.y[i] %in% c('%', 'percent')){
                        if (!is.na(sd.y[i])){
                            # good
                        }else{
                            sd.val = sd13 * 100 # 7JUL20
                            sd.y[i] = sd.val
                        }
                        units.y[i] = '%'
                        
                    }else if (units.y[i] %in% c('mol/kg', 'mmol/g')){
                        new.val = val.y[i] * .001
                        if (!is.na(sd.y[i])){
                            sd.val = sd.y[i] * .001
                            sd.y[i] = sd.val
                        }else{
                            sd.val = new.val * sd16 # 7JUL20
                            sd.y[i] = sd.val
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mol/g'
                        
                    }else if (units.y[i] %in% c('mmol/kg', 'micromol/g')){
                        new.val = val.y[i] * .000001
                        if (!is.na(sd.y[i])){
                            sd.val = sd.y[i] * .000001
                            sd.y[i] = sd.val
                        }else{
                            sd.val = new.val * sd16 # 7JUL20
                            sd.y[i] = sd.val
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mol/g'
                        
                    }else if (units.y[i] == 'mol/m2'){
                        new.val = val.y[i] * areas[i]
                        if (!is.na(sd.y[i])){
                            sd.val = sd.y[i] * areas[i]
                            sd.y[i] = sd.val
                        }else{
                            sd.val = new.val * sd16 # 7JUL20
                            sd.y[i] = sd.val
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mol/g'
                        
                    }else if (units.y[i] == 'mol/g'){
                        if (!is.na(sd.y[i])){
                            # good
                        }else{
                            sd.val = val.y[i] * sd16 # 7JUL20
                            sd.y[i] = sd.val
                        }
                        
                    }else if (units.y[i] == 'fraction'){
                        if (!is.na(sd.y[i])){
                            # good
                        }else{
                            sd.val = sd13 # 7JUL20
                            sd.y[i] = sd.val
                        }
                        
                    }else if (units.y[i] == 'micromol/m2'){ # 19JUL20
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        
                        new.val = val.y[i] * (1 / 1000000) * areas[i] * min.val
                        if (!is.na(sd.y[i])){
                            new.sd = sd.y[i] * (1 / 1000000) * areas[i] * min.val
                            sd.y[i] = new.sd
                        }else{
                            new.sd = new.val * sd16
                            sd.y[i] = new.sd
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mol/L'
                    }
                    
                }else if (axis.y[i] == 'aqueous'){
                    molar.mass = mass(s_element(sorbents[i]))
                    if (units.y[i] == 'ppb'){
                        new.val = (val.y[i] * 0.001) / (1000 * molar.mass)
                        if (!is.na(sd.y[i])){
                            new.sd = (sd.y[i] * 0.001) / (1000 * molar.mass)
                            sd.y[i] = new.sd
                        }else{
                            new.sd = new.val * sd16 # 7JUL20
                            sd.y[i] = new.sd
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mol/L'
                        
                    }else if (units.y[i] == 'ppm'){
                        if (sorbents[i] == 'Zn(+2)'){
                            molar.mass = 65.38
                        }else{
                            molar.mass = mass(s_element(sorbents[i]))
                        }
                        new.val = (val.y[i]) / (1000 * molar.mass)
                        if (!is.na(sd.y[i])){
                            new.sd = (sd.y[i]) / (1000 * molar.mass)
                            sd.y[i] = new.sd
                        }else{
                            new.sd = new.val * sd16 # 7JUL20
                            sd.y[i] = new.sd
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mol/L'
                        
                    }else if (units.y[i] == 'mol/L'){
                        if (!is.na(sd.y[i])){
                            # good
                        }else{
                            new.sd = val.y[i] * sd16 # 7JUL20
                            sd.y[i] = new.sd
                        }
                    }
                    
                }else if (axis.y[i] == 'Kd' | axis.y[i] == 'log(Kd)'){
                    if (units.y[i] %in% c('mL/g', 'ml/g')){
                        units.y[i] = 'mL/g'
                    }
                    
                }else if (axis.y[i] == 'Rd' | axis.y[i] == 'log(Rd)'){
                    if (units.y[i] == 'L/kg'){
                        units.y[i] = 'mL/g'
                        
                    }else if (units.y[i] == 'm3/kg'){
                        new.val = val.y[i] * 1000
                        if (!is.na(sd.y[i])){
                            sd.val = sd.y[i] * 1000
                            sd.y[i] = sd.val
                        }
                        val.y[i] = new.val
                        units.y[i] = 'mL/g'
                    }
                }
            }
        }
        
        dataset$Y_axis = axis.y # update values
        dataset$Y_val = val.y # update values
        dataset$Y_SD = sd.y  # update values
        dataset$Y_units = units.y # update values
        
        #############################################################
        ### HANDLE X-AXIS CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        #############################################################
        axis.x = dataset$X_axis # read values
        val.x = dataset$X_val # read values
        sd.x = dataset$X_SD # read values
        units.x = dataset$X_units # read values
        
        minerals = dataset$Mineral # read values
        val.minerals = dataset$Mineral_val # read values
        sd.minerals = dataset$Mineral_SD # read values
        units.minerals = dataset$Mineral_units # read values
        
        # clean X data values
        for (i in c(1:p)){
            if (units.x[i] != '' & !is.na(units.x[i])){
                if (axis.x[i] == 'pH'){ # transfer pH before calculations, moved from top of X cleaning
                    axis.x[i] = NA
                    val.x[i] = NA
                    sd.x[i] = NA
                    units.x[i] = NA
                }else if (axis.x[i] == 'aqueous'){
                    if (units.x[i] == 'log(mol/L)'){
                        new.val = 10^(val.x[i])
                        if (!is.na(sd.x[i])){
                            new.sd = 10^(sd.x[i])
                            sd.x[i] = new.sd
                        }else{
                            new.sd = new.val * sd17 # 7JUL20
                            sd.x[i] = new.sd
                        }
                        val.x[i] = new.val
                        axis.x[i] = 'aqueous'
                        units.x[i] = 'mol/L'
                        
                    }else if (units.x[i] == 'mol/L'){
                        if (!is.na(sd.x[i])){
                            # good
                        }else{
                            new.sd = val.x[i] * sd16 # 7JUL20
                            sd.x[i] = new.sd
                        }
                    }else if (units.x[i] == 'mmol/L'){
                        new.val = val.x[i] * (1 / 1000)
                        if (!is.na(sd.x[i])){
                            new.sd = sd.x[i] * (1 / 1000)
                            sd.x[i] = new.sd
                        }else{
                            new.sd = new.val * sd16 # 7JUL20
                            sd.x[i] = new.sd
                        }
                        val.x[i] = new.val
                        axis.x[i] = 'aqueous'
                        units.x[i] = 'mol/L'
                    }
                    
                }else if (axis.x[i] == 'log(aq)'){
                    new.val = 10^(val.x[i])
                    if (!is.na(sd.x[i])){
                        new.sd = 10^(sd.x[i])
                        sd.x[i] = new.sd
                    }else{
                        new.sd = new.val * sd17 # 7JUL20
                        sd.x[i] = new.sd
                    }
                    val.x[i] = new.val
                    axis.x[i] = 'aqueous'
                    units.x[i] = 'mol/L'
                    
                }else if (axis.x[i] == 'aq'){
                    if (units.x[i] == 'mmol/L'){
                        new.val = val.x[i] / 1000
                        if (!is.na(sd.x[i])){
                            new.sd = sd.x[i] / 1000
                            sd.x[i] = new.sd
                        }else{
                            new.sd = new.val * sd16 # 7JUL20
                            sd.x[i] = new.sd
                        }
                        val.x[i] = new.val
                        axis.x[i] = 'aqueous'
                        units.x[i] = 'mol/L'
                    }
                    
                }else if (axis.x[i] == 'loading'){
                    if (units.x[i] == 'mol/kg'){
                        new.val = val.x[i] / 1000
                        if (!is.na(sd.x[i])){
                            new.sd = sd.x[i] / 1000
                            sd.x[i] = new.sd
                        }else{
                            new.sd = new.val * sd16 # 7JUL20
                            sd.x[i] = new.sd
                        }
                        val.x[i] = new.val
                        axis.x[i] = 'sorbed'
                        units.x[i] = 'mol/g'
                    }
                    
                }else if (axis.x[i] %in% c('bentonite', 'mineral', 'solid', 'mineral_val')){ # 19JUL20
                    molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                    new.val = val.x[i] / molar.mass
                    if (!is.na(sd.x[i])){
                        new.sd = sd.x[i] / molar.mass
                        sd.minerals[i] = new.sd
                    }else{
                        new.sd = new.val * sd5 # 7JUL20
                        sd.minerals[i] = new.sd
                    }
                    val.minerals[i] = new.val
                    units.minerals[i] = 'mol/L'
                    axis.x[i] = NA
                    val.x[i] = NA
                    sd.x[i] = NA
                    units.x[i] = NA
                    
                }else if (axis.x[i] == 'sorbed'){
                    if (units.x[i] == 'mol/kg'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        new.val = (val.x[i] / 1000) * min.val
                        if (!is.na(sd.x[i])){
                            new.sd = (sd.x[i] / 1000) * min.val
                            sd.x[i] = new.sd
                        }else{
                            new.sd = new.val * sd16 # 7JUL20 # 19JUL20
                            sd.x[i] = new.sd
                        }
                        val.x[i] = new.val
                        units.x[i] = 'mol/L'
                    }
                    
                    ### SURFACE CHARGE CONVERSIONS AND TRANSFERS ###
                }else if (axis.x[i] == 'total H(+)' | axis.x[i] == 'H(+)'){ # 8JUL20
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                    }else{
                        molar.mass = mass(s_element(minerals[i]))
                    }
                    min.val = val.minerals[i] * molar.mass
                    
                    if (units.x[i] == 'mol/L'){
                        new.val = (val.x[i] - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                        if (!is.na(sd.x[i])){
                            new.sd = (sd.x[i] - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 16JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                        axis.x[i] = 'charge'
                        
                    }else if (units.x[i] == 'mmol/L'){ # 19JUL20
                        new.val = ((val.x[i]/1000) - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                        if (!is.na(sd.x[i])){
                            new.sd = ((sd.x[i]/1000) - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 16JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                        axis.x[i] = 'charge'
                        
                    }
                    
                }else if (axis.x[i] == 'total OH(-)' | axis.x[i] == 'OH(-)'){
                    if (units.x[i] == 'mol/L'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        
                        new.val = (-1 * val.x[i] - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                        if (!is.na(sd.x[i])){
                            new.sd = (-1 * sd.x[i] - (10^(-pHs[i])) + (10^(-1*(14-pHs[i])))) * (1/min.val) * (1/areas[i]) * (1000000/10000) * 96485
                            sd.charges[i] = new.sd
                        }else{
                            new.sd = new.val * sd11 # 16JUL20
                            sd.charges[i] = new.sd
                        }
                        val.charges[i] = new.val
                        units.charges[i] = 'microC/cm2'
                        axis.x[i] = 'charge'
                    }
                    
                }else if (units.x[i] == 'microC/cm2'){
                    if (!is.na(sd.x[i])){
                        sd.charges[i] = sd.x[i]
                    }else{
                        new.sd = val.x[i] * sd11 # 7JUL20
                        sd.charges[i] = new.sd
                    }
                    val.charges[i] = new.val
                    units.charges[i] = 'microC/cm2'
                    axis.x[i] = 'charge'
                }
                
                # transfer alkalinity (HCO3-) and elytes to elytes columns
                if (axis.x[i] %in% c('pH', 'aqueous', 'sorbed', 'charge', '') | is.na(axis.x[i])){
                    # do nothing, these are good
                    
                }else if (axis.x[i] == 'alkalinity'){
                    for (n in c(1:6)){
                        elyte = paste('Electrolyte', n, sep = '')
                        elyte.val = paste('Electrolyte', n, '_val', sep = '')
                        elyte.sd = paste('Electrolyte', n, '_SD', sep = '')
                        elyte.unit = paste('Electrolyte', n, '_units', sep = '')
                        elytes = as.vector(dataset[,elyte])
                        val.elytes = as.vector(dataset[,elyte.val])
                        sd.elytes = as.vector(dataset[,elyte.sd])
                        units.elytes = as.vector(dataset[,elyte.unit])
                        
                        if (is.na(elytes[i]) | elytes[i] == ''){
                            elytes[i] = 'HCO3(-)' # transfer elyte name
                            val.elytes[i] = val.x[i] / 1000 # transfer elyte val
                            if (!is.na(sd.elytes[i])){
                                sd.elytes[i] = sd.x[i] / 1000 # transfer elyte sd
                            }else{
                                sd.elytes[i] = val.elytes[i] * sd2 # 7JUL20
                            }
                            units.elytes[i] = 'mol/L' # transfer elyte units
                            axis.x[i] = NA
                            val.x[i] = NA
                            sd.x[i] = NA
                            units.x[i] = NA
                            dataset[,elyte] = elytes
                            dataset[,elyte.val] = val.elytes
                            dataset[,elyte.sd] = sd.elytes
                            dataset[,elyte.unit] = units.elytes
                            break
                            
                        }
                    }
                }else{
                    for (n in c(1:6)){
                        elyte = paste('Electrolyte', n, sep = '')
                        elyte.val = paste('Electrolyte', n, '_val', sep = '')
                        elyte.sd = paste('Electrolyte', n, '_SD', sep = '')
                        elyte.unit = paste('Electrolyte', n, '_units', sep = '')
                        elytes = as.vector(dataset[,elyte])
                        val.elytes = as.vector(dataset[,elyte.val])
                        sd.elytes = as.vector(dataset[,elyte.sd])
                        units.elytes = as.vector(dataset[,elyte.unit])
                        
                        if (is.na(elytes[i]) | elytes[i] == ''){
                            if (units.x[i] == 'mg/L'){ # 7JUL20
                                if (axis.x[i] %in% mineral.ref$minerals){ 
                                    molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == axis.x[i], 'masses'])
                                    new.val = val.x[i] / (1000 * molar.mass)
                                }
                                elytes[i] = axis.x[i] # transfer elyte name
                                val.elytes[i] = new.val # transfer elyte val # 19JUL20
                                if (!is.na(sd.elytes[i])){
                                    sd.elytes[i] = sd.x[i] / (1000 * molar.mass) # transfer elyte sd
                                }else{
                                    sd.elytes[i] = val.elytes[i] * sd2 # 7JUL20
                                }
                                units.elytes[i] = 'mol/L' # transfer elyte units
                            }else{
                                elytes[i] = axis.x[i] # transfer elyte name
                                val.elytes[i] = val.x[i] # transfer elyte val
                                if (!is.na(sd.elytes[i])){
                                    sd.elytes[i] = sd.x[i] # transfer elyte sd
                                }else{
                                    sd.elytes[i] = val.elytes[i] * sd2 # 7JUL20
                                }
                                units.elytes[i] = units.x[i] # transfer elyte units
                            }
                            axis.x[i] = NA
                            val.x[i] = NA
                            sd.x[i] = NA
                            units.x[i] = NA
                            dataset[,elyte] = elytes
                            dataset[,elyte.val] = val.elytes
                            dataset[,elyte.sd] = sd.elytes
                            dataset[,elyte.unit] = units.elytes
                            break
                            
                        }
                    }
                }
            }
        }
        
        dataset$X_axis = axis.x # update values
        dataset$X_val = val.x # update values
        dataset$X_SD = sd.x # update values
        dataset$X_units = units.x # update values
        
        dataset$Mineral_val = val.minerals # update values
        dataset$Mineral_SD = sd.minerals # update values
        dataset$Mineral_units = units.minerals # update values
        
        dataset$SurfCharge_val = val.charges # update values
        dataset$SurfCharge_SD = sd.charges # update values
        dataset$SurfCharge_units = units.charges # update values
        
        dataset$pH = pHs
        dataset$pH_SD = sd.pHs
        
        #####################################################################
        ### HANDLE SORBED/AQUEOUS CONVERSIONS, TRANSFERS, AND ESTIMATIONS ###
        #####################################################################
        sorbents = dataset$Sorbent # read values
        val.sorbents = dataset$Sorbent_val # read values
        sd.sorbents = dataset$Sorbent_SD # read values
        units.sorbents = dataset$Sorbent_units # read values
        
        val.sorbed = dataset$Sorbed_val # read values
        sd.sorbed = dataset$Sorbed_SD # read values
        units.sorbed = dataset$Sorbed_units # read values
        
        val.aqueous = dataset$Aq_val # read values
        sd.aqueous = dataset$Aq_SD # read values
        units.aqueous = dataset$Aq_units # read values
        
        axis.x = dataset$X_axis # read values
        axis.x = replace(axis.x, c(axis.x == ''), NA)
        val.x = dataset$X_val # read values
        sd.x = dataset$X_SD # read values
        units.x = dataset$X_units # read values
        
        axis.y = dataset$Y_axis # read values
        axis.y = replace(axis.y, c(axis.y == ''), NA)
        val.y = dataset$Y_val # read values
        sd.y = dataset$Y_SD # read values
        units.y = dataset$Y_units # read values
        
        minerals = dataset$Mineral # read values
        val.minerals = dataset$Mineral_val # read values
        
        for (i in c(1:p)){
            if (!is.na(axis.x[i]) & !is.na(axis.y[i])){ # cases where neither X or Y is NA
                
                if (axis.x[i] == 'aqueous' & axis.y[i] == 'sorbed'){
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                    }else{
                        molar.mass = mass(s_element(minerals[i]))
                    }
                    min.val = val.minerals[i] * molar.mass
                    
                    if (units.y[i] == 'mol/mol'){ # 19JUL20
                        val.sorbed[i] = val.y[i] * val.minerals[i]
                        sd.sorbed[i] = sd.y[i] * val.minerals[i]
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.x[i]
                        sd.aqueous[i] = sd.x[i]
                        units.aqueous[i] = 'mol/L'
                        val.sorbents[i] = val.sorbed[i] + val.aqueous[i]
                        sd.sorbents[i] = val.sorbents[i] * sd16
                        units.sorbents[i] = 'mol/L'
                        
                    }else{
                        val.sorbed[i] = val.y[i] * min.val # val.minerals[i] removed on 16JUL20
                        sd.sorbed[i] = sd.y[i] * min.val # val.minerals[i] removed on 16JUL20
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.x[i]
                        sd.aqueous[i] = sd.x[i]
                        units.aqueous[i] = 'mol/L'
                        val.sorbents[i] = val.sorbed[i] + val.aqueous[i]
                        sd.sorbents[i] = val.sorbents[i] * sd16 # 7JUL20
                        units.sorbents[i] = 'mol/L'
                    }
                    
                }else if (axis.x[i] == 'sorbed' & axis.y[i] == 'Kd'){ # 9JUL20
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                    }else{
                        molar.mass = mass(s_element(minerals[i]))
                    }
                    min.val = val.minerals[i] * molar.mass
                    val.aqueous[i] = val.x[i] / val.y[i]
                    sd.aqueous[i] = val.aqueous[i] - (val.x[i] / (val.y[i]) * (1 + sd14))
                    units.aqueous[i] = 'mol/L'
                    val.sorbed[i] = val.x[i]
                    sd.sorbed[i] = sd.x[i]
                    units.sorbed[i] = 'mol/L'
                    val.sorbents[i] = val.aqueous[i] + val.sorbed[i]
                    sd.sorbents[i] = val.sorbents[i] * sd10
                    units.sorbents[i] = 'mol/L'
                    
                }else if (axis.x[i] == 'aqueous' & axis.y[i] == 'Rd'){
                    if (units.x[i] == 'mol/L'){
                        val.sorbed[i] = val.sorbents[i] - val.x[i]
                        sd.sorbed[i] = sd.x[i]
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.sorbents[i] - val.sorbed[i]
                        sd.aqueous[i] = sd.x[i]
                        units.aqueous[i] = 'mol/L'
                    }
                }else if (axis.x[i] == 'aqueous' & axis.y[i] == 'log(Kd)'){
                    # handle pbb99_f6b case here
                }
                
            }else if (is.na(axis.x[i]) & !is.na(axis.y[i])){ # cases where only X is NA
                if (axis.y[i] == 'aqueous'){
                    val.sorbed[i] = val.sorbents[i] - val.y[i]
                    sd.sorbed[i] = sd.y[i]
                    units.sorbed[i] = 'mol/L'
                    val.aqueous[i] = val.y[i]
                    sd.aqueous[i] = sd.y[i]
                    units.aqueous[i] = 'mol/L'
                    
                }else if (axis.y[i] == 'sorbed'){
                    if (units.y[i] == 'mol/g'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        
                        val.sorbed[i] = val.y[i] * min.val # val.minerals[i] removed on 16JUL20
                        sd.sorbed[i] = sd.y[i] * min.val # val.minerals[i] removed on 16JUL20
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.sorbents[i] - val.sorbed[i]
                        sd.aqueous[i] = sd.y[i] * val.minerals[i]
                        units.aqueous[i] = 'mol/L'
                        
                    }else if (units.y[i] == '%'){
                        val.sorbed[i] = val.sorbents[i] * (val.y[i] / 100)
                        sd.sorbed[i] = val.sorbents[i] * sd13 # 7JUL20
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.sorbents[i] - val.sorbed[i]
                        sd.aqueous[i] = val.sorbents[i] * sd13 # 7JUL20
                        units.aqueous[i] = 'mol/L'
                        
                    }else if (units.y[i] == 'fraction'){
                        val.sorbed[i] = val.sorbents[i] * val.y[i]
                        sd.sorbed[i] = val.sorbents[i] * sd13 # 7JUL20
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.sorbents[i] - val.sorbed[i]
                        sd.aqueous[i] = val.sorbents[i] * sd13 # 7JUL20
                        units.aqueous[i] = 'mol/L'
                        
                    }else if (units.y[i] == 'mol/L'){ # 19JUL20
                        val.sorbed[i] = val.y[i]
                        sd.sorbed[i] = sd.y[i] # 7JUL20
                        units.sorbed[i] = 'mol/L'
                        val.aqueous[i] = val.sorbents[i] - val.sorbed[i]
                        sd.aqueous[i] = sd.y[i] # 7JUL20
                        units.aqueous[i] = 'mol/L'
                    }
                    
                }else if (axis.y[i] == 'Rd' | axis.y[i] == 'Kd'){
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                    }else{
                        molar.mass = mass(s_element(minerals[i]))
                    }
                    min.val = val.minerals[i] * molar.mass
                    val.aqueous[i] = val.sorbents[i] / ((val.y[i] * (1/1000) * min.val + 1))
                    sd.aqueous[i] = (((val.aqueous[i] - (val.sorbents[i] / ((val.y[i]*(1+sd14)) / 1000 * min.val + 1))) / val.aqueous[i])^2 + sd10^2)^0.5 * val.aqueous[i]
                    units.aqueous[i] = 'mol/L'
                    val.sorbed[i] = val.sorbents[i] - val.aqueous[i]
                    sd.sorbed[i] = (sd.sorbents[i]^2 + sd.aqueous[i]^2)^(0.5)
                    units.sorbed[i] = 'mol/L'
                    
                }else if (axis.y[i] == 'log(Rd)' | axis.y[i] == 'log(Kd)'){
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                    }else{
                        molar.mass = mass(s_element(minerals[i]))
                    }
                    min.val = val.minerals[i] * molar.mass
                    val.aqueous[i] = val.sorbents[i] / ((val.y[i] * (1/1000) * min.val + 1))
                    sd.aqueous[i] = (((val.aqueous[i]- (val.sorbents[i] / (10^(log(val.y[i])+sd15) / 1000 * min.val + 1))) / val.aqueous[i])^2 + sd10^2)^0.5 * val.aqueous[i]
                    units.aqueous[i] = 'mol/L'
                    val.sorbed[i] = val.sorbents[i] - val.aqueous[i]
                    sd.sorbed[i] = (sd.sorbents[i]^2 + sd.aqueous[i]^2)^(0.5)
                    units.sorbed[i] = 'mol/L'
                }else if (axis.y[i] == 'Ka'){ # 8JUL20
                    if (units.y[i] == 'mL/m2'){
                        if (minerals[i] %in% mineral.ref$minerals){
                            molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        }else{
                            molar.mass = mass(s_element(minerals[i]))
                        }
                        min.val = val.minerals[i] * molar.mass
                        val.aqueous[i] = (val.sorbents[i] / ((val.y[i] * areas[i]) * (1 / 1000) * min.val + 1))
                        sd.aqueous[i] = val.aqueous[i] - (val.sorbents[i] / ((val.y[i] * (1 + sd14) * areas[i]) * (1 / 1000) * min.val + 1))
                        units.aqueous[i] = 'mol/L'
                        val.sorbed[i] = val.sorbents[i] - val.aqueous[i]
                        sd.sorbed[i] = (sd.sorbents[i]^2 + sd.aqueous[i]^2)^(0.5)
                        units.sorbed[i] = 'mol/L'
                    }
                }
                
            }else if (!is.na(axis.x[i]) & is.na(axis.y[i])){ # cases where only Y is NA
                # handle only x-axis case
            }
        }
        
        dataset$Sorbed_val = val.sorbed # update values
        dataset$Sorbed_SD = sd.sorbed # update values
        dataset$Sorbed_units = units.sorbed # update values
        
        dataset$Aq_val = val.aqueous # update values
        dataset$Aq_SD = sd.aqueous # update values
        dataset$Aq_units = units.aqueous # update values
        
        dataset$Sorbent_val = val.sorbents # update values
        dataset$Sorbent_SD = sd.sorbents # update values
        dataset$Sorbent_units = units.sorbents
        
        ########################################
        ### CONVERT MINERAL INFO BACK TO G/L ###
        ########################################
        val.minerals = dataset$Mineral_val # read values
        sd.minerals = dataset$Mineral_SD # read values
        units.minerals = dataset$Mineral_units # read values
        
        for (i in c(1:p)){
            if (units.minerals[i] != '' & !is.na(units.minerals[i])){
                if (units.minerals[i] == 'mol/L'){
                    if (minerals[i] %in% mineral.ref$minerals){
                        molar.mass = as.numeric(mineral.ref[mineral.ref$minerals == minerals[i], 'masses'])
                        new.val = val.minerals[i] * molar.mass
                        if (!is.na(sd.minerals[i])){
                            new.sd = sd.minerals[i] * molar.mass
                            sd.minerals[i] = new.sd
                        }
                        val.minerals[i] = new.val
                        units.minerals[i] = 'g/L'
                    }
                }
            }
        }
        
        dataset$Mineral_val = val.minerals # update values
        dataset$Mineral_SD = sd.minerals # update values
        dataset$Mineral_units = units.minerals # update values
        
        ######################################
        ### REMOVE X/Y COLUMNS FROM OUTPUT ###
        ######################################
        col.keep = c('X_axis', 'X_val', 'X_SD', 'X_units', 'Y_axis', 'Y_val', 'Y_SD', 'Y_units')
        dataset = dataset[,!(colnames(dataset) %in% col.keep)]
        
        ####################################
        ### CONVERT MISSING VALUES TO NA ###
        ####################################
        for (i in c(1:ncol(dataset))){
            dirty = as.vector(dataset[,i])
            # clean = na_if(dirty, 0)
            cleanest = na_if(dirty, '')
            dataset[,i] = cleanest
        }
        
        ###############################################
        ### ROUND ALL NUMERIC COLUMNS TO 6 SIG FIGS ###
        ###############################################
        for (i in c(1:ncol(dataset))){
            if (is.double(dataset[,i])){
                values = as.vector(dataset[,i])
                dataset[,i] = signif(values, 6)
            }
        }
        
        ##############
        ### OUTPUT ###
        ##############
        dataset
    })
    
    output$sc.cleaned = DT::renderDataTable({ # output datatable to UI
        sc.cleaned = sc.dataset()
        sc.cleaned
    })
    
    output$downloadData = downloadHandler(
        filename = function() {
            paste('sc.dataset', '.csv', sep = '')
        },
        content = function(file) {
            write.csv(sc.dataset(), file, row.names = F) # user must import in Excel, not just open
    })
})
