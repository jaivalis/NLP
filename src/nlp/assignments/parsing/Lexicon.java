package nlp.assignments.parsing;

import nlp.ling.Tree;
import nlp.util.Counter;
import nlp.util.CounterMap;

import java.util.List;
import java.util.Set;

/**
 * Simple default implementation of a lexicon, which scores word, tag pairs with a smoothed estimate of
 * P(tag|word)/P(tag).
 */
class Lexicon {
    CounterMap<String, String> wordToTagCounters = new CounterMap<String, String>();
    double totalTokens = 0.0;
    double totalWordTypes = 0.0;
    Counter<String> tagCounter = new Counter<String>();
    Counter<String> wordCounter = new Counter<String>();
    Counter<String> typeTagCounter = new Counter<String>();

    public Set<String> getAllTags() {
        return tagCounter.keySet();
    }

    public boolean isKnown(String word) {
        return wordCounter.keySet().contains(word);
    }

    public double scoreTagging(String word, String tag) {
        double p_tag = tagCounter.getCount(tag) / totalTokens;
        double c_word = wordCounter.getCount(word);
        double c_tag_and_word = wordToTagCounters.getCount(word, tag);
        if (c_word < 10) { // rare or unknown
            c_word += 1.0;
            // closed word classes will receive very low c_tag_and_word
            c_tag_and_word += typeTagCounter.getCount(tag) / totalWordTypes;
        }
        // add plus one smoothing
        double p_word = (1.0 + c_word) / (totalTokens + 1.0);
        double p_tag_given_word = c_tag_and_word / c_word;
        return p_tag_given_word / p_tag * p_word;
    }

    public Lexicon(List<Tree<String>> trainTrees) {
        for (Tree<String> trainTree : trainTrees) {
            List<String> words = trainTree.getYield();
            List<String> tags = trainTree.getPreTerminalYield();
            for (int position = 0; position < words.size(); position++) {
                String word = words.get(position);
                String tag = tags.get(position);
                tallyTagging(word, tag);
            }
        }
    }

    private void tallyTagging(String word, String tag) {
        if (!isKnown(word)) {
            totalWordTypes += 1.0;
            typeTagCounter.incrementCount(tag, 1.0);
        }
        totalTokens += 1.0;
        tagCounter.incrementCount(tag, 1.0);
        wordCounter.incrementCount(word, 1.0);
        wordToTagCounters.incrementCount(word, tag, 1.0);
    }
}