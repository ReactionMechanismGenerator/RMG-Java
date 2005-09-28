GCJ = jc#/usr/local/bin/gcj #java compiler
OBJECTS = GasConstant.o InvalidUnitException.o ParameterInfor.o Pressure.o Temperature.o ChemParser.o InvalidGraphFormatException.o InvalidReactionAdjListFormatException.o InvalidReactionFormatException.o InvalidThermoFormatException.o InvalidUnionFormatException.o NotInDictionaryException.o ReplaceFunctionalGroupException.o Atom.o AtomOcuppiedException.o BondDictionary.o Bond.o ChemElementDictionary.o ChemElement.o ChemGraph.o ChemNodeElement.o DuplicatedElementException.o EmptyAtomException.o FailGenerateThermoDataException.o FGAtom.o FGElementDictionary.o FGElement.o ForbiddenStructureException.o FreeElectronDictionary.o FreeElectron.o FunctionalGroupCollection.o FunctionalGroupNotFoundException.o FunctionalGroup.o GATP.o GeneralGAPP.o GroupNotFoundException.o InvalidAtomListException.o InvalidBondException.o InvalidBondMutationException.o InvalidCenterTypeException.o InvalidChangedOrderException.o InvalidChemElementException.o InvalidChemGraphException.o InvalidChemNodeElementException.o InvalidDelocalizationPathException.o InvalidFreeElectronException.o InvalidFunctionalGroupException.o InvalidNodeElementException.o InvalidReferenceThermoGAValueException.o InvalidSpeciesException.o InvalidThermoCenterException.o LennardJones.o Matchable.o MultipleGroupFoundException.o NASAThermoData.o NullSymbolException.o ReplaceThermoGAValueException.o SiteNotInSpeciesException.o SpeciesDictionary.o Species.o TemperatureOutOfRangeException.o ThermoData.o ThermoGAGroupLibrary.o ThermoGAValue.o ThreeFrequencyModel.o UnknownSpeciesNameException.o UnknownSymbolException.o Arc.o chemUtil_pkgClass.o DummyLeaf.o FailToIdentifyFGAtomException.o GraphComponent.o Graph.o HierarchyTreeNodeException.o HierarchyTreeNode.o HierarchyTree.o InconnectedGraphException.o InvalidCentralIDException.o InvalidChildException.o InvalidConnectivityException.o InvalidCycleDetectionException.o InvalidGraphComponentException.o InvalidHierarchyRelationException.o InvalidMatchedSiteException.o InvalidNeighborException.o InvalidNeighborSizeException.o InvalidNodeIDException.o InvalidTreeNodeLevelException.o MatchedSite.o Node.o NotInGraphException.o NotPartitionedException.o NullGraphComponentException.o NullGraphException.o OverConnectedException.o PositionOccupiedException.o ReplaceCentralNodeException.o TreeNode.o Tree.o Curve.o EmptyQueueException.o FullQueueException.o InvalidDoubleFormatException.o InvalidIntegerFormatException.o InvalidUncertaintyException.o InvalidUncertaintyTypeException.o MathTool.o Queue.o UncertainDouble.o Action.o ArrheniusEPKinetics.o ArrheniusKinetics.o DuplicatedIsomersException.o FailToGenerateReverseReactionException.o InvalidActionException.o InvalidDistanceException.o InvalidKineticsException.o InvalidKineticsFormatException.o InvalidKineticsKeyException.o InvalidKineticsTemplateException.o InvalidKineticsTypeException.o InvalidPDepNetReactionException.o InvalidPDepNetworkTypeException.o InvalidProductException.o InvalidProductNumberException.o InvalidReactantException.o InvalidReactantNumberException.o InvalidReactantTreeException.o InvalidReactionDirectionException.o InvalidReactionTemplateDirectionException.o InvalidStructureException.o InvalidStructureTemplateException.o InvalidTemplateReactionException.o Kinetics.o KineticsTemplateLibrary.o KineticsTemplate.o LibraryReactionGenerator.o LibraryReaction.o MatchedLeafNotFoundException.o MultipleReverseReactionException.o NegativeAException.o NegativeRateException.o PDepNetReaction.o PDepNetwork.o PDepPathReaction.o PDepReactionTemplate.o PDepWell.o RateConstantNotFoundException.o RateConstant.o ReactionAdjList.o ReactionGenerator.o ReactionLibrary.o Reaction.o ReactionTemplateLibrary.o ReactionTemplate.o ReplaceKineticsTemplateException.o SpeciesTemplateTree.o Structure.o StructureTemplate.o StructureTree.o TemplateReactionGenerator.o TemplateReaction.o ThirdBodyReaction.o TROEReaction.o AbstractReactionModel.o AdiabaticTM.o Chemkin.o ConstantPM.o ConstantTM.o ConversionTT.o CoreEdgeReactionModel.o Core.o CurvedPM.o CurvedTM.o DAESolver.o DynamicSimulatorException.o DynamicSimulator.o Edge.o FinishController.o GenericTM.o InitialStatus.o InvalidBeginStatusException.o InvalidChemkinParameterException.o InvalidConversionException.o InvalidCoreEdgeRelationException.o InvalidNextCandidateSpeciesException.o InvalidPressureModelException.o InvalidReactedReactionException.o InvalidReactedSpeciesException.o InvalidReactionModelEnlargerException.o InvalidReactionModelTypeException.o InvalidReactionSetException.o InvalidReactionSystemUpdateException.o InvalidReactionTimeException.o InvalidReactionTimeUnitException.o InvalidSpeciesStatusException.o InvalidSymbolException.o InvalidSystemCompositionException.o InvalidTemperatureModelException.o InvalidUnreactedReactionException.o JDASPK.o JDASSL.o Milestone.o NegativeConcentrationException.o ODEReaction.o ODESolver.o PresentStatus.o PressureModel.o PrimaryReactionLibrary.o RateBasedPDepRME.o RateBasedPDepVT.o RateBasedRME.o RateBasedVT.o Reactant.o ReactionModelEnlarger.o ReactionModelGenerator.o ReactionModel.o ReactionNotAddedException.o ReactionSystem.o ReactionTime.o ReactionTimeTT.o SASolver.o SpeciesConversion.o SpeciesStatus.o SystemSnapshot.o TemperatureModel.o TerminationTester.o UnknownReactedSpeciesException.o ValidityTester.o RMG.o Maths.o CholeskyDecomposition.o EigenvalueDecomposition.o LUDecomposition.o Matrix.o QRDecomposition.o SingularValueDecomposition.o 

%.o: jing/param/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/chemParser/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/chem/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/chemUtil/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/mathTool/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/param/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/rxn/%.java
	$(GCJ) -c -g -classpath . $<
%.o:jing/rxnSys/%.java
	$(GCJ) -c -g -classpath . $<
%.o:Jama/%.java
	$(GCJ) -c -g $<
%.o:Jama/util/%.java
	$(GCJ) -c -g $<

.PHONY: all
all: $(OBJECTS) RMG.exec

JDASPK.o:jing/rxnSys/JDASPK.java
	$(GCJ) -c -g -fjni -classpath . $<
RMG.o:RMG.java
	$(GCJ) -c -g $<
RMG.exec:$(OBJECTS)
	$(GCJ) -o RMG.exec --main=RMG $(OBJECTS)

.PHONY: clean
clean: 
	rm *.o
