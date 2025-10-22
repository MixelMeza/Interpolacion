
from flask import Flask, render_template, request
from sympy import symbols, expand, simplify

app = Flask(__name__)

def lagrange_interpolacion(x, xs, ys):
    n = len(xs)
    resultado = 0
    for i in range(n):
        termino = ys[i]
        for j in range(n):
            if i != j:
                termino *= (x - xs[j]) / (xs[i] - xs[j])
        resultado += termino
    return resultado

# Polinomio simbólico
def lagrange_polinomio(xs, ys):
    x = symbols('x')
    n = len(xs)
    polinomio = 0
    for i in range(n):
        Li = 1
        for j in range(n):
            if i != j:
                Li *= (x - xs[j])/(xs[i] - xs[j])
        polinomio += ys[i]*Li
    return simplify(expand(polinomio))

@app.route('/', methods=['GET', 'POST'])
def index():
    resultado = None
    polinomio = None
    puntos = []
    x_eval = ''
    if request.method == 'POST':
        try:
            xs = [float(x) for x in request.form.getlist('x')]
            ys = [float(y) for y in request.form.getlist('y')]
            x_eval = float(request.form['x_eval'])
            resultado = lagrange_interpolacion(x_eval, xs, ys)
            polinomio = lagrange_polinomio(xs, ys)
            puntos = list(zip(xs, ys))
        except Exception:
            resultado = 'Error en los datos ingresados.'
    return render_template('lagrange_form.html', resultado=resultado, polinomio=polinomio, puntos=puntos, x_eval=x_eval)

if __name__ == '__main__':
    app.run(debug=True)
