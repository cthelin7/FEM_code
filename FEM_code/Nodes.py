import numpy

nodes_list = {}


class Node:
    """

    """

    def __init__(self, node_id):
        """

        """

        self.node_id = node_id
        self.location = []
        self.value = None
        # nodes_list.append(self)
        nodes_list[self.node_id] = self

    def add_location(self, location):
        self.location = location
