import numpy

nodes_list = []


class Node:
    """

    """

    def __init__(self, id):
        """

        """

        self.node_id = id
        self.location = []
        self.value = None
        nodes_list.append(self)

    def add_location(self, location):
        self.location = location
