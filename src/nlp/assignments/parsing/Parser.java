package nlp.assignments.parsing;

import nlp.ling.Tree;

import java.util.List;

/**
 * Parsers are required to map sentences to trees.  How a parser is constructed and trained is not specified.
 */
interface Parser {
    Tree<String> getBestParse(List<String> sentence);

    double getLogScore(Tree<String> tree);
}