import sys
import random


def sample_multiple_queries(community_file, query_file, test_file):
    total_communities = []
    community_f = open(community_file, 'r')
    for line in community_f:
        total_communities.append(line)
    community_f.close()

    community_orders = set()
    while len(community_orders) < 1000:
        candidate = random.randint(0, len(total_communities) - 1)
        community_orders.add(candidate)

    query_f = open(query_file, 'w')
    test_f = open(test_file, 'w')

    counter = 0
    for community_order in community_orders:
        community = total_communities[community_order]
        members = community.strip().split()

        query_number = random.randint(1, 16)
        while query_number > len(members):
            query_number = random.randint(1, 16)

        query_vertices = set()
        while len(query_vertices) < query_number:
            query_vertices.add(members[random.randint(0, len(members) - 1)])

        query_string = ''
        for query_vertex in query_vertices:
            query_string = query_string + query_vertex + ' '

        query_string = query_string[:-1] + '\n'

        query_f.write(query_string)
        test_f.write(community)
        counter += 1
        print 'multiple community generated: ' + str(counter)

    query_f.close()
    test_f.close()


if __name__ == '__main__':

    community_file = sys.argv[1]
    query_file = sys.argv[2]
    test_file = sys.argv[3]
    sample_multiple_queries(community_file, query_file, test_file)
