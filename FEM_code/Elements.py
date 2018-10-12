import numpy

elements_list = []


class Element:
    """

    """
    def __init__(self, id):
        """

        """

        self.element_id = id
        self.nodes = []

        elements_list.append(self)

    def add_nodes(self, node_ids):
        self.nodes.append(node_ids)
