import ast
import operator as op
from dataclasses import dataclass



@dataclass
class Workspace:
    """
    Workspace for mathematical expressions.
    """

    functions: dict

    OPERATORS = {
        ast.Add:  op.add,
        ast.Sub:  op.sub,
        ast.Mult: op.mul,
        ast.Div:  op.truediv,
        ast.Pow:  op.pow,
        ast.USub: op.neg
        }


    def eval(self, node, names={}, functions={}):
        """
        Evaluate (previously parsed) expression based on dictionary of *names* and *functions*.
        """
        N, F = names, functions

        # number
        if isinstance(node, ast.Num):
            return node.n

        # constant
        elif isinstance(node, ast.Constant):
            return node.value

        # tuples
        elif isinstance(node, ast.Tuple):
            return tuple(self.eval(x, N, F) for x in node.elts)

        # names
        elif isinstance(node, ast.Name):
            if not node.id in N:
                raise(NameError(f"name '{node.id}' is not defined"))
            return N[node.id]

        # binary operator
        elif isinstance(node, ast.BinOp):
            left, right = self.eval(node.left, N, F), self.eval(node.right, N, F)
            return self.OPERATORS[type(node.op)](left, right)

        # unary operator
        elif isinstance(node, ast.UnaryOp):
            return self.OPERATORS[type(node.op)](self.eval(node.operand, N, F))

        # function evaluation
        elif isinstance(node, ast.Call):
            fid = node.func.id
            if fid in self.functions:
                func = self.functions[fid]
            elif fid in F:
                func = F[fid]
            else:
                raise(NameError(f"unkown or unsupported function '{fid}'"))

            args = [self.eval(arg, N, F) for arg in node.args]
            return func(*args)

        # type error
        else:
            raise TypeError(node)


    def parse(self, expression):
        """
        Parse expression into AST node.
        """
        return ast.parse(expression, mode='eval')


    def dependencies(self, expression):
        """
        List of dependencies for computing variable expression.
        """
        root = self.parse(expression)
        required = lambda node: isinstance(node, ast.Name) and not node.id in self.functions
        return sorted([node.id for node in ast.walk(root) if required(node)])
