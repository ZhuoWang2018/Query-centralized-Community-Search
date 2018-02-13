import os
import sys


def compute_f1_score(test_c, result_c):
    test_set = set(test_c)
    result_set = set(result_c)
    intersect = len(result_set.intersection(test_set))
    prec = (intersect + 0.0) / len(result_set)
    recall = (intersect + 0.0) / len(test_set)
    f1 = 2 * prec * recall / (prec + recall)
    return f1, prec, recall


def handle_multi_communities(test_c, result_f):
    max_f1 = -1.0
    max_prec = -1.0
    max_recall = -1.0
    for line in result_f:
        if line.startswith('#'):
            break
        else:
            result_c = line.strip().split()
            f1, prec, recall = compute_f1_score(test_c, result_c)
            if max_f1 < f1:
                max_f1 = f1
                max_prec = prec
                max_recall = recall
    return max_f1, max_prec, max_recall


def measure_f1_score(test_file, result_file, query_file):
    test_f = open(test_file, 'r')
    query_f = open(query_file, 'r')
    result_f = open(result_file, 'r')

    f1 = 0.0
    prec = 0.0
    recall = 0.0

    counter = 0
    for line in result_f:
        test_c = test_f.readline().strip().split()
	query_c = query_f.readline().strip().split()
	skip = False
	if len(query_c) == 1:
	    skip = True

        if line.startswith('#'):
            delta_f1, delta_prec, delta_recall = handle_multi_communities(test_c, result_f)
	    if skip == False:
                f1 += delta_f1
                prec += delta_prec
                recall += delta_recall
        else:
            result_c = line.strip().split()
            delta_f1, delta_prec, delta_recall = compute_f1_score(test_c, result_c)
	    if skip == False:
                f1 += delta_f1
                prec += delta_prec
                recall += delta_recall
	if skip == False:
            counter += 1
    test_f.close()
    query_f.close()
    result_f.close()
    print 'counter:' + str(counter)
    return f1/counter, prec/counter, recall/counter

if __name__ == '__main__':


    if len(sys.argv) < 4:
        print 'test file!'
        print 'query file!'
        print 'result file!'
	os._exit(0)

    test_file = sys.argv[1]

    query_file = sys.argv[2]

    result_file = sys.argv[3]



    print test_file
    print result_file

    f1, prec, recall = measure_f1_score(test_file, result_file, query_file)
    print 'precision: ' + str(prec)
    print 'recall: ' + str(recall)
    print 'f1: ' + str(f1)
