package nlp.assignments.parsing;

import java.util.*;

import nlp.ling.Tree;
import nlp.util.CounterMap;


class BaselineCkyParser implements Parser {

	  
    CounterMap<List<String>, Tree<String>> knownParses;
    CounterMap<Integer, String> spanToCategories;
    
    TreeAnnotator annotator;
    
    Lexicon lexicon;
    Grammar grammar;
    
    //UnaryClosure unaryClosure;
    
    static class Chart {
    	/*
    	 * TODO
    	 *  This class (and enclosed EdgeInfo) needs to be changed to keep track of unary rules (and unary chains) used
    	 */
    	
        static class EdgeInfo {
        	double score; 
        	BinaryRule rule = null;
        	int mid = -1;
        	EdgeInfo(double score) {
        		this.score = score;
        	}
			@Override
			public String toString() {
				if (rule == null) {
					return Double.toString(score);
				} else {
					return score + ": " + "[ rule = " + rule + ", mid = " + mid + "]";
				}
				
			}
        }

    	Map<Integer,Map<Integer,Map<String, EdgeInfo>>> chart = new HashMap<Integer, Map<Integer,Map<String,EdgeInfo>>>();
    	Chart(int seqLength) {
    		for (int i = 0; i < seqLength; i++) {
    			chart.put(i, new HashMap<Integer, Map<String,EdgeInfo>>());
    			for (int j = i + 1; j <= seqLength; j++) {
    				chart.get(i).put(j, new HashMap<String, EdgeInfo>());
    			}
    		}
    	}
    	
    	void set(int i, int j, String label, double score) {
    		chart.get(i).get(j).put(label, new EdgeInfo(score));
    	}
    	double get(int i, int j, String label) {
    		Map<String,EdgeInfo> edgeScores = chart.get(i).get(j);
    		if (!edgeScores.containsKey(label)) {
    			return 0;
    		} else {
    			return edgeScores.get(label).score;
    		}
    	}
    	
    	Set<String> getAllCandidateLabels(int i, int j) {
    		return chart.get(i).get(j).keySet();
    	}
    	
    	String getBestLabel(int i, int j) {
    		Map<String,EdgeInfo> edgeScores = chart.get(i).get(j);
    		double bestScore = Double.NEGATIVE_INFINITY;
    		String optLabel = null;
    		for (String label : edgeScores.keySet()) {
    			if (bestScore < edgeScores.get(label).score) {
    				optLabel = label;
    				bestScore = edgeScores.get(label).score;
    			}
    		}
    		return optLabel;
    	}
    	
    	void setBackPointer(int i, int j, String label, BinaryRule rule, int midPoint) {
    		EdgeInfo edgeInfo = chart.get(i).get(j).get(label);
    		edgeInfo.rule = rule;
    		edgeInfo.mid = midPoint;
    	}

    	
    	
    	int getMidPoint(int i, int j, String label) {
    		return chart.get(i).get(j).get(label).mid;
    	}

    	BinaryRule getRule(int i, int j, String label) {
    		return chart.get(i).get(j).get(label).rule;
    	}

		@Override
		public String toString() {
			return chart.toString();
		}

    }
    
    void traverseBackPointersHelper(List<String> sent, Chart chart, int i, int j, Tree<String> currTree) {
    	String parent = currTree.getLabel();
    	
    	/** 
    	 * TODO 
    	 * This method needs to be updated to keep print out unary rules used 
    	 */
    	
		// binary rules
    	if (j - i > 1) {

    		BinaryRule rule = chart.getRule(i, j, parent);
       		int mid = chart.getMidPoint(i, j, parent);
       		List<Tree<String>> children = new ArrayList<Tree<String>>(2);
    		
    		Tree<String> t1 = new Tree<String>(rule.getLeftChild()); 
    		traverseBackPointersHelper(sent, chart, i, mid, t1);
    		children.add(t1);
    		
    		Tree<String> t2 = new Tree<String>(rule.getRightChild());
    		traverseBackPointersHelper(sent, chart, mid, j, t2);
    		children.add(t2);
    		
    		currTree.setChildren(children);
    		
    	// preterminal production
    	} else {
    		assert j - i == 1;
    		
    		Tree<String> termProd = new Tree<String>(sent.get(i));
    		currTree.setChildren(Collections.singletonList(termProd));
    	}
    	
    }

	// traverse back pointers and create a tree
    Tree<String> traverseBackPointers(List<String> sentence, Chart chart) {

    	Tree<String> annotatedBestParse;
    	if (!chart.getAllCandidateLabels(0, sentence.size()).contains("ROOT")) {
    		// this line is here only to make sure that a baseline without binary rules can output something 
    		annotatedBestParse = new Tree<String>(chart.getBestLabel(0, sentence.size()));
    	} else {
    		// in reality we always want to start with the ROOT symbol of the grammar
    		annotatedBestParse = new Tree<String>("ROOT");
    	}
        traverseBackPointersHelper(sentence, chart, 0, sentence.size(), annotatedBestParse);
        return annotatedBestParse;
    	
    }
    
    
    public Tree<String> getBestParse(List<String> sentence) {
    	
      /*
       * TODO
       * This method needs to be extended to support unary rules 
       * The UnaryClosure class should simplify this task substantially 
       */

      // note that chart.get(i, j, c) translates to chart[i][j][c] we used in the slides	
    	
      Chart chart = new Chart(sentence.size());
      
      // preterminal rules
      for (int k = 0; k < sentence.size(); k++) {
    	  for (String preterm : lexicon.getAllTags()) {
    		  chart.set(k, k + 1, preterm, lexicon.scoreTagging(sentence.get(k), preterm));
    		  
    	  }
      }

      // CKY for binary trees
      for (int max = 2; max <= sentence.size(); max++) {
    	  for (int min = max - 2; min >= 0; min--) {
    		  // first without unary roles
    		  for (String parent : grammar.states) {
    			  double bestScore = Double.NEGATIVE_INFINITY;
    			  int optMid =-1;
    			  BinaryRule optRule = null;
    			  // parent -> c1 c2
    			  for (BinaryRule rule : grammar.getBinaryRulesByParent(parent)) {
    				  for (int mid = min + 1; mid <  max; mid++) {
    					  double score1 = chart.get(min, mid, rule.getLeftChild());
    					  double score2 = chart.get(mid, max, rule.getRightChild());
    					  double currScore = score1 * score2 * rule.getScore();
    					  if (currScore > bestScore) {
    						  bestScore = currScore;
    						  optMid = mid;
    						  optRule = rule;;
    					  }
    				  }
    			  }
    			  if (bestScore != Double.NEGATIVE_INFINITY) { 
    				  chart.set(min, max, parent, bestScore);
					  chart.setBackPointer(min, max, parent, optRule, optMid);
				  }
    		  }
    		  
    	  }
      }
       
      // use back pointers to create a tree
      Tree<String> annotatedBestParse = traverseBackPointers(sentence, chart);

      return annotator.unAnnotateTree(annotatedBestParse);
    }

   

    public BaselineCkyParser(List<Tree<String>> trainTrees, TreeAnnotator annotator) {

      this.annotator = annotator;

    	
      System.out.print("Annotating / binarizing training trees ... ");
      List<Tree<String>> annotatedTrainTrees = annotateTrees(trainTrees);

      System.out.println("done.");

      System.out.print("Building grammar ... ");
      grammar = new Grammar(annotatedTrainTrees);
      System.out.print("Building grammar ... ");
      lexicon = new Lexicon(annotatedTrainTrees);
      
      //this.lexicon.scoreTagging(word, tag)
      System.out.println("done. (" + grammar.getStates().size() + " states)");
      
      // use the unary closure to support unary rules in the CKY algorithm
      //unaryClosure = new UnaryClosure(grammar);
    
    }

    private List<Tree<String>> annotateTrees(List<Tree<String>> trees) {
      List<Tree<String>> annotatedTrees = new ArrayList<Tree<String>>();
      for (Tree<String> tree : trees) {
    	  annotatedTrees.add(annotator.annotateTree(tree));
      }
      return annotatedTrees;
    }

    /**
     * Recursive function that populates a given collection with the Unary rules from a given tree.
     * @param tree The tree to be processed.
     * @param ret A collection to be populated.
     */
    private void extractUnaryRules(Tree<String> tree, Collection<UnaryRule> ret) {
        if (tree.getChildren().size() == 0) { return; } // termination condition
        else if (tree.getChildren().size() == 1) {      // 1 child > add to ret
            if (!tree.getChildren().get(0).isLeaf()){
                String parent = tree.getLabel();
                String child = tree.getChildren().get(0).getLabel();

                ret.add(new UnaryRule(parent, child));
            }
            extractUnaryRules(tree.getChildren().get(0), ret);
        }
        else if  (tree.getChildren().size() == 2) {
            extractUnaryRules(tree.getChildren().get(0), ret);
            extractUnaryRules(tree.getChildren().get(1), ret);
        }
    }

    /**
     * Recursive function that populates a given collection with the pre-terminal Unary rules from a given tree.
     * @param tree The tree to be processed.
     * @param ret A collection to be populated.
     */
    private void extractPreterminalRules(Tree<String> tree, Collection<UnaryRule> ret) {
        if (tree.getChildren().size() == 0) { return; } // termination condition
        else if (tree.getChildren().size() == 1) {      // 1 child > add to ret
            if (tree.isPreTerminal()) {
                String parent = tree.getLabel();
                String child = tree.getChildren().get(0).getLabel();

                ret.add(new UnaryRule(parent, child));
            }
            extractPreterminalRules(tree.getChildren().get(0), ret);
        }
        else if  (tree.getChildren().size() == 2) {
            extractPreterminalRules(tree.getChildren().get(0), ret);
            extractPreterminalRules(tree.getChildren().get(1), ret);
        }
    }

    /**
     * Recursive function that populates a given collection with the binary rules from a given tree.
     * @param tree The tree to be processed.
     * @param ret A collection to be populated.
     */
    private void extractBinaryRules(Tree<String> tree, Collection<BinaryRule> ret) {
        if (tree.getChildren().size() == 0) { return; }
        else if (tree.getChildren().size() == 1) {
            extractBinaryRules(tree.getChildren().get(0), ret);
        }
        else if  (tree.getChildren().size() == 2) {
            String parent = tree.getLabel();
            String left = tree.getChildren().get(0).getLabel();
            String right = tree.getChildren().get(1).getLabel();
            ret.add(new BinaryRule(parent, left, right));

            extractBinaryRules(tree.getChildren().get(0), ret);
            extractBinaryRules(tree.getChildren().get(1), ret);
        }
    }
    
    @Override
    public double getLogScore(Tree<String> tree) {
        double score = 0.0;
        Tree<String> annotatedTree = annotator.annotateTree(tree);
        /*
         * Add method which uses the annotatedTree (not the 'tree') and compute the log probability of the tree
         */
        Collection<UnaryRule> preterminalRules = new ArrayList<UnaryRule>();
        extractPreterminalRules(annotatedTree, preterminalRules);

        Collection<UnaryRule> unaryRules = new ArrayList<UnaryRule>();
        extractUnaryRules(annotatedTree, unaryRules);

        Collection<BinaryRule> binaryRules = new ArrayList<BinaryRule>();
        extractBinaryRules(annotatedTree, binaryRules);


        for (UnaryRule pr : preterminalRules) {
            double sc = this.lexicon.scoreTagging(pr.child, pr.parent);
            System.out.println(pr);
            score += Math.log(sc);
        }
        for (UnaryRule ur : unaryRules) {
            double sc = this.lexicon.scoreTagging(ur.child, ur.parent);
            //System.out.println(ur);
            //System.out.println(sc);
            score += Math.log(sc);
        }
        for (BinaryRule br : binaryRules) {
            int index = this.grammar.getBinaryRules().indexOf(br);
            double ruleLogP = Math.log(0);
            if (index != -1) {
                ruleLogP = Math.log(this.grammar.getBinaryRules().get(index).getScore());
            }
            score += ruleLogP;
        }
        return score;
    }


}