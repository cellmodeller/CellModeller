from libsbml import * 
import math
import new
import urllib2


# indicates whether we should perform SBML document consistency check.
# supressed by default because the MSR SBML library does not produce valid SBML. 
supressConsistencyCheck = True


# get the array index of given species ID:
def getSpecIndex(specIdToArrayIndexDict, specId):
    if not specIdToArrayIndexDict.has_key(specId):
        specIdToArrayIndexDict[specId] = len(specIdToArrayIndexDict.keys())    
    return specIdToArrayIndexDict[specId]       


# get python variable name for the ODE of given species ID: 
def specOdeNameFromId(specId):
    return "d_" + specId


# checks the given SBML document for consistency and raises
# an exception if internal inconsistency is detected.
def checkSBMLConsistency(document):
    if supressConsistencyCheck:
        return
    
    numFailures = document.checkInternalConsistency()
    if numFailures > 0:
        failureMsg = ""
        for failureNum in range(numFailures):
            print "Failure " + str(failureNum) + ": " + document.getError(failureNum).getShortMessage() + "\n"
        raise Exception(failureMsg)
            
            
# load SBML model from given file.
def SBMLModelFromSBMLFile(sbmlFile):    
    reader = SBMLReader()    
    document = reader.readSBML(sbmlFile)
    if document.getNumErrors()>0:
        print "Errors in reading SBML file..."
    checkSBMLConsistency(document)
    model = document.getModel()
    if not model:
        print "No model!"
    return model


# load SBML model form given string.
def SBMLModelFromSBMLStr(sbmlStr):
    
    # would be nice if SBMLReader.readSBMLFromStr worked - but it does not, 
    # there seems to be some SWIG fuck-up. so we are forced to first write
    # to temporary file:     
    file = open("tmp.sbml", 'w')
    file.write(sbmlStr) 
    file.close()    
    
    # then read from file:
    reader = SBMLReader()    
    document = reader.readSBML("tmp.sbml")
    checkSBMLConsistency(document)    
    model = document.getModel()
    
    return model


# helper function for composing a python string from a binary AST math node.
# mutually recursive with pythonMathFromASTNode.
def splitBinary(astNode, operand, kineticLaw, model):
    str1 = pythonMathFromASTNode(astNode.getLeftChild(), kineticLaw, model)
    str2 = pythonMathFromASTNode(astNode.getRightChild(), kineticLaw, model)
    return "(" + str1 + " " + operand + " " + str2 + ")"


# constructs a python string representing the math expresson in an AST node.
# currently only covers basic cases, and not all the ASTNode types; an 
# exception is thrown if an unsupported type is encountered.
# more types can easily be added.
def pythonMathFromASTNode(astNode, kineticLaw, model):
    if (astNode.getType() == AST_PLUS):
        return splitBinary(astNode, "+", kineticLaw, model)
        
    elif (astNode.getType() == AST_MINUS):
        return splitBinary(astNode, "-", kineticLaw, model)

    elif (astNode.getType() == AST_DIVIDE):
        return splitBinary(astNode, "/", kineticLaw, model)
                    
    elif (astNode.getType() == AST_TIMES):
        return splitBinary(astNode, "*", kineticLaw, model)             
        
    elif (astNode.getType() == AST_FUNCTION_POWER):
        return splitBinary(astNode, "**", kineticLaw, model)             

    elif astNode.isNumber():
        if astNode.isReal():
            return str(astNode.getReal())
        elif astNode.isInteger():
            return str(astNode.getInteger())    
        
    elif astNode.isName():
        # check if this name is a defined parameter. can't see any way of distinguishing between
        # parameters and species identifiers, other than checking for definition this way?
        nodeName = astNode.getName()
        parameterLocal = kineticLaw.getParameter(nodeName)
        parameterGlobal = model.getParameter(nodeName)               
        if parameterLocal != None:
            val = str(parameterLocal.getValue())            
            return val
        elif parameterGlobal != None:
            val = str(parameterGlobal.getValue())
            return val
        
        # if not parameter, then check if this is a species id:    
        elif model.getSpecies(nodeName) != None:           
            return nodeName
        
        # something weird -- return 1 and issue warning:
        else:
            print "WARNING: unknown identifier " + nodeName + ", defaulting to 1. Globals: " + str(model.getNumParameters()) + ". Locals: " + str(kineticLaw.getNumParameters())
            return "1"
        
    else:    
        raise Exception("Un-supported AST node type: " + str(astNode.isName()) + ", node: " + str(astNode.getType()))
        

# constructs a python program string from given SBML model.
def pythonStrFromSBMLModel(model):

    # initialise dictionaries and python code string:    
    specIdToArrayIndexDict = {}      
    specIdToExpDict = {}
    
    # define indentation for python code:
    indent = "    "
    indent2 = indent + indent   
    
    # iterate over reactions:
    for reactionNum in range(model.getNumReactions()):
  
        reaction = model.getReaction(reactionNum)
        
        # compute rate expression. 
        # UPDATE: first tried to use the formulaToString libSBML function,
        # which should have worked in most cases. but it doesn't substitute
        # in defined parameters. hence we use our own ASTNode parsing function.    
        #mathAstNode = reaction.getKineticLaw().getMath()            
        #pythonRateExp = "(" + formulaToString(mathAstNode) + ")"
        #pythonRateExp = reaction.getKineticLaw().getFormula()        
        pythonRateExp = pythonMathFromASTNode(reaction.getKineticLaw().getMath(), reaction.getKineticLaw(), model)         
                                    
        """
        # iterate over reactants to create mass-action rate expression:
        # UPDATE: don't assume mass-action, just read the rate expression literally.
        for reactantNum in range(reaction.getNumReactants()):
            
            reactantSpecRef = reaction.getReactant(reactantNum)
            stoich = reactantSpecRef.getStoichiometry()
            
            # for some reason stoichiometry is sometimes NaN:
            if math.isnan(stoich):
                stoich = 1
                    
            reactantSpecId = reactantSpecRef.getSpecies()
                    
            if stoich == 1:
                pythonRateExp += "*" + reactantSpecId
            else: 
                pythonRateExp += "*" + reactantSpecId + "^" + str(stoich)
        """
                               
        # iterate over reactants and update ODE for each:      
        for reactantNum in range(reaction.getNumReactants()):
            reactantSpecRef = reaction.getReactant(reactantNum)
            stoich = reactantSpecRef.getStoichiometry()
            
            # for some reason stoichiometry is sometimes NaN:
            if math.isnan(stoich):
                stoich = 1
            
            reactantSpecId = reactantSpecRef.getSpecies()        
            
            # add key to hash table if it doesn't exist:       
            if not specIdToExpDict.has_key(reactantSpecId):       
                specIdToExpDict[reactantSpecId] = ""
                   
            if stoich == 1:
                specIdToExpDict[reactantSpecId] += " - " + pythonRateExp
            else:
                specIdToExpDict[reactantSpecId] += " - " + str(stoich) + "*" + pythonRateExp
            
        # do the same for products:    
        for productNum in range(reaction.getNumProducts()):
            productSpecRef = reaction.getProduct(productNum)
            stoich = productSpecRef.getStoichiometry()
            
            # for some reason stoichiometry is sometimes NaN:
            if math.isnan(stoich):
                stoich = 1
            
            productSpecId = productSpecRef.getSpecies()        
                   
            # add key to hash table if it doesn't exist:       
            if not specIdToExpDict.has_key(productSpecId):       
                specIdToExpDict[productSpecId] = ""       
                   
            if stoich == 1:
                specIdToExpDict[productSpecId] += " + " + pythonRateExp 
            else:
                specIdToExpDict[productSpecId] += " + " + str(stoich) + "*" + pythonRateExp
                  
    # create python variable definitions, python updated rate definitions,
    # and a list of species IDs in the relevant order.
    pythonVarDefs = ""
    pythonRateUpdateExps = ""        
    specIdLst = [0]*len(specIdToExpDict)
    for specId in specIdToExpDict.keys():
        # get index into species array:
        index = getSpecIndex(specIdToArrayIndexDict, specId)
        
        # create python variable definition:
        pythonVarDefs += indent + indent + specId + " = cell.species[" + str(index) + "]\n"

        # create python rate update expression:
        exp = specIdToExpDict[specId]
        varname = specOdeNameFromId(specId)
        pythonRateUpdateExps += indent2 + varname + " = " + exp + "\n"
        
        # insert spec id into correct position in list:
        specIdLst[index] = specId #.insert(index, specId)        
            
    # create a python list with updated rate expressions,
    # and python lists with with species IDs and names in the relevant order.      
    pythonRateUpdateLst = "" 
    pythonSpecIDLst = ""
    pythonSpecNameLst = ""
    for specId in specIdLst:
        
        # get species name from model:               
        specName = model.getSpecies(specId).getName()
        if specName == None:
            specName = "Undefined"            
        
        if pythonRateUpdateLst == "":
            pythonRateUpdateLst += "[" + specOdeNameFromId(specId)
            pythonSpecIDLst += "[\"" + specId + "\""
            pythonSpecNameLst += "[\"" + specName + "\""
        else:
            pythonRateUpdateLst += ", " + specOdeNameFromId(specId)
            pythonSpecIDLst += ", \"" + specId + "\""
            pythonSpecNameLst += ", \"" + specName + "\""
            
    pythonRateUpdateLst += "]"
    pythonSpecIDLst += "]"
    pythonSpecNameLst += "]"          
    
    # build initial population list:
    initPopLst = [0]*model.getNumSpecies()
    for specId in specIdLst:        
        sbmlSpec = model.getSpecies(specId)
        initAmount = sbmlSpec.getInitialAmount()
        initConc = sbmlSpec.getInitialConcentration()
        init = 0.0
        if not math.isnan(initAmount) and initAmount > 0:
            init = initAmount
        elif not math.isnan(initConc) and initConc > 0:
            init = initConc
            
        index = getSpecIndex(specIdToArrayIndexDict, specId)    
        initPopLst[index] = init
        
    # build initial population python string:
    pythonInitPopLst = ""
    for init in initPopLst:
        if pythonInitPopLst == "":
            pythonInitPopLst += "[" + str(init)
        else:
            pythonInitPopLst += "," + str(init)
            
    pythonInitPopLst += "]"    
    
    # compose the complete python program string.    
    pythonStr = ""
    pythonStr += "def getRates(cells):\n"
    pythonStr += indent + "df = []\n"
    pythonStr += indent + "for cell in cells:\n" 
    pythonStr += pythonVarDefs + "\n"
    pythonStr += pythonRateUpdateExps + "\n"
    pythonStr += indent2 + "df.append(" + pythonRateUpdateLst + ")\n" 
    pythonStr += indent + "return df"        
    pythonStr += "\n\n"
    pythonStr += "specIdLst = " + pythonSpecIDLst + "\n\n"
    pythonStr += "specNameLst = " + pythonSpecNameLst + "\n\n"
    pythonStr += "specInitPopLst = " + pythonInitPopLst + "\n\n"
        
    return pythonStr

    
# create python module from code dynamically.    
def pythonModuleFromPythonStr(pythonStr):
    module = new.module("sbmlPythonEncoding")
    exec pythonStr in module.__dict__
    return module    


# get text from given channel.    
def TextFromChannel(channelName):
    f = urllib2.urlopen('http://async-message-passer.appspot.com/?content_only=1&channel_name=' + channelName)
    return f.read()    


def pythonModuleFromFile(fileName):
    reader = SBMLReader()    
    document = reader.readSBML(fileName)
    checkSBMLConsistency(document)    
    model = document.getModel()   
    pstr = pythonStrFromSBMLModel(model)
    print "Got python str: " + pstr
    pythonModule = pythonModuleFromPythonStr(pstr)
    return pythonModule


# get python module from given channel.
def pythonModuleFromChannel(channelName):
    sbmlStr = TextFromChannel(channelName)

    #print "Got channel text:\n" + sbmlStr          
    
    # 1) an apparent segmentation issue in libsbml means that we cannot
    # invoke SBMLModelFromSBMLStr() defined above. great. so we inline it below:
    # 2) would be nice if SBMLReader.readSBMLFromStr worked - but it does not, 
    # there seems to be some SWIG fuck-up. so we are forced to first write
    # to temporary file:
    file = open("tmp.sbml", 'w') 
    file.write(sbmlStr) 
    file.close()    
    reader = SBMLReader()    
    document = reader.readSBML("tmp.sbml")
    checkSBMLConsistency(document)    
    model = document.getModel()   
        
                
    pythonStr = pythonStrFromSBMLModel(model)
    print "Got python str: " + pythonStr
    pythonModule = pythonModuleFromPythonStr(pythonStr)
    return pythonModule


# example usage:    
#sbmlFile = "/Users/michaeldavidpedersen/tmp/testsbml.sbml"  
#print sbmlToPython(sbmlFile)
#module = pythonModuleFromChannel("DEFAULT_CHANNEL")
