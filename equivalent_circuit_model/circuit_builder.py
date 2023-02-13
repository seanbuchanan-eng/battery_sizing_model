class Circuit:
	"""
	This class represents an equivalent circuit battery model.
	It is a mutable object.

	Attributes
	----------

	Methods
	-------
	""" 

	def __init__(self):
		""" 
		Parameters
		----------
		"""
		self._circuit_elements = []
		# do

	def __str__(self):
		""" 
		Parameters
		----------
		"""
		node = 0
		print_string = ""
		for element in self._circuit_elements:
			print_string += f"Node {node}: {element}\n"
			node += 1

		return print_string

	def print_equation(self):
		"""
		Print the equations for this circuit as pertaining to 
		Kirkoff's laws.

		Parameters
		----------
		Nothing

		Returns
		-------
		str
			a string representation of the equation representing this circuit.


		"""

	def add_open_circuit_voltage(self, connection):
		"""
		Add an open circuit voltage element to the circuit.

		This element is added to the node directly preceding it.

		Parameters
		----------
		connection : str
			Either "series" or "parallel". Indicates the type of connection to the 
			previously added element. 

		Returns
		-------
		Nothing
		"""
		print("ocv")

	def add_resistor(self, name, connection, resistance=0):
		"""
		Add a resistor element to the circuit.

		This element is added to the node directly preceding it.

		Parameters
		----------
		name : str
			Name of the resistor in the circuit. Must be unique.
		connection : str
			Either "series" or "parallel". Indicates the type of connection to the 
			previously added element.
		resistance: float
			Resistance of the resistor (default is 0).

		Raises
		-------
		ValueError
			If connection is not "series" or "parallel"
		"""
		resistor = Resistor(name, resistance)

		if (connection == "series"):
			circuit_node = CircuitNode("series", [resistor])
			self._circuit_elements.append(circuit_node)
		elif (connection == "parallel"):
			if (len(self._circuit_elements) < 1):
				raise ValueError('First circuit element must be a series element.')
			last_element = self._circuit_elements[-1].get_elements()
			del self._circuit_elements[-1]
			self._circuit_elements.append(CircuitNode("parallel", last_element.append(resistor)))
		else:
			raise ValueError('connection must be "series" or "parallel"')

	def add_capacitor(self, connection):
		"""
		Add a capacitor element to the circuit.

		This element is added to the node directly preceding it.

		Parameters
		----------
		connection : str
			Either "series" or "parallel". Indicates the type of connection to the 
			previously added element. 

		Returns
		-------
		Nothing


		"""
		print("capacitor")
	def add_inductor(self, connection):
		"""
		Add an inductor element to the circuit.

		This element is added to the node directly preceding it.

		Parameters
		----------
		connection : str
			Either "series" or "parallel". Indicates the type of connection to the 
			previously added element. 

		Returns
		-------
		Nothing


		"""
		print("inductor")

class CircuitNode:
	"""
	This class represents a node in a circuit.

	A CircuitNode defines whether a connection between circuit elements
	is series or parallel. If the connection is parallel it contains all
	of the elements within that parallel section.
	"""

	def __init__(self, connection, elements):
		"""
		Parameters
		----------
		connection : str
			The type of connection of the node. Either "series" or "parallel".
		elements : list
			The elements wihin the node section. A length of 1 for series and > 1 for parallel nodes.
		"""

		self.elements = elements
		self.connection = connection

	def __str__(self):
		return ','.join(str(element) for element in self.elements)

	def is_parallel(self):
		"""
		Returns
		-------
		bool
			true if a parellel node else false
		"""
		return self.connection == "parallel"

	def get_elements(self):
		"""
		Returns
		-------
		list
			List of circuit elements in the node
		"""
		return self.elements

class Resistor:
	"""
	Represents a resistor in an electrical circuit.
	"""

	def __init__(self, name, resistance):
		"""
		Parameters
		----------
		name : str
			Name of the resistor.
		resitance : float
			Resistance of the resitor in Ohms.
		"""
		self.name = name
		self.resistance = resistance

	def __str__(self):
		return f"Resistor {self.name} {self.resistance} Ohms"

	def get_name(self):
		"""
		Returns
		-------
		str
			Name of the resistor.
		"""
		return self.name

	def get_resistance(self):
		"""
		Returns
		-------
		float
			Resitance value of the resistor.
		"""
		return self.resistance

class Capacitor:
	"""
	Represents a capacitor in an electrical circuit.
	"""

	def __init__(self, capacitance):
		"""
		Parameters
		----------
		capacitance : float
			Capacitance of the capacitor in farads.
		"""
		self.capacitance = capacitance

	def get_capacitance(self):
		"""
		Returns
		-------
		float
			Capacitance value of the capacitor.
		"""
		return self.capacitance

class OpenCircuitVoltage:
	"""
	Represents an open circuit voltage element in an electrical circuit.
	"""

	def __init__(self):
		self.voltage = None

	def get_voltage(self, SOC):
		"""
		Parameters
		----------
		SOC : float (0, 1]
			The current state of charge of the battery.

		Returns
		-------
		float
			Open circuit voltage of the circuit at SOC.
		"""
		print("unimplimented")

class Inductor:
	"""
	Represents an inductor in an electrical circuit.
	"""

	def __init__(self, inductance):
		"""
		Parameters
		----------
		inductance : float
			Inductance value of the inductor in Henry.
		"""
		self.inductance = inductance

	def get_inductance(self):
		"""
		Returns
		-------
		float
			Inductance value of the inductor in Henry.
		"""
		return self.inductance