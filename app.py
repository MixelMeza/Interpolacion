
from flask import Flask, render_template, request
from sympy import symbols, expand, simplify

app = Flask(__name__)


def lagrange_interpolacion(x, xs, ys):
    n = len(xs)
    resultado = 0
    pasos = []
    for i in range(n):
        termino = ys[i]
        paso = f"L_{i}(x) = {ys[i]}"
        for j in range(n):
            if i != j:
                termino *= (x - xs[j]) / (xs[i] - xs[j])
                paso += f" * (({x} - {xs[j]}) / ({xs[i]} - {xs[j]}))"
        resultado += termino
        paso += f" = {termino}"
        pasos.append(paso)
    pasos.append(f"Suma total: {resultado}")
    return resultado, pasos

def lagrange_polinomio(xs, ys):
    x = symbols('x')
    n = len(xs)
    polinomio = 0
    pasos = []
    for i in range(n):
        Li = 1
        paso = f"L_{i}(x) = "
        for j in range(n):
            if i != j:
                Li *= (x - xs[j])/(xs[i] - xs[j])
                paso += f"((x - {xs[j]})/({xs[i]} - {xs[j]}))"
        polinomio += ys[i]*Li
        paso = f"{ys[i]}*" + paso if paso != f"L_{i}(x) = " else f"{ys[i]}"
        pasos.append(paso)
    pasos.append(f"Polinomio expandido: {simplify(expand(polinomio))}")
    return simplify(expand(polinomio)), pasos

# Newton: diferencias divididas

def tabla_diferencias_divididas(xs, ys):
    n = len(xs)
    # Para visualización (shift visual hacia abajo)
    tabla = [[None for _ in range(n+2)] for _ in range(n)]
    # Para coeficientes (sin shift, solo fila 0)
    dd = [ys[:]]
    debug_msgs = []
    try:
        # Llenar i, xi, f(xi)
        for i in range(n):
            tabla[i][0] = i
            tabla[i][1] = xs[i]
            tabla[i][2] = ys[i]
        # Llenar diferencias divididas para visualización (shift visual hacia abajo)
        for orden in range(1, n):
            for fila in range(orden, n):
                arriba = tabla[fila][orden+1-1]
                abajo = tabla[fila-1][orden+1-1]
                divisor = xs[fila] - xs[fila-orden]
                debug_msgs.append(f"orden={orden}, fila={fila}, arriba={arriba}, abajo={abajo}, divisor={divisor}")
                if divisor == 0:
                    valor = None
                else:
                    valor = (arriba - abajo) / divisor
                tabla[fila][orden+2-1] = valor
        # Llenar diferencias divididas para coeficientes (sin shift, solo fila 0)
        for j in range(1, n):
            col = []
            for i in range(n-j):
                divisor = xs[i+j] - xs[i]
                if divisor == 0:
                    col.append(None)
                else:
                    col.append((dd[j-1][i+1] - dd[j-1][i]) / divisor)
            dd.append(col)
    except Exception as e:
        raise Exception(f"Error en tabla_diferencias_divididas: {e}. Debug: {debug_msgs}")
    return tabla, dd

def tabla_html(tabla_dd):
    # tabla_dd must be a tuple (tabla_visual, dd)
    if not isinstance(tabla_dd, tuple) or len(tabla_dd) != 2:
        return "<div>Error: tabla interna no disponible para visualización.</div>"
    tabla_vis, dd = tabla_dd
    n = len(tabla_vis)
    max_orden = len(dd) - 1
    headers = ['i', 'xᵢ', 'f(xᵢ)'] + [f"Orden {j}" for j in range(1, max_orden+1)]
    html = "<b>Tabla de diferencias divididas:</b>"
    html += "<table border='1' style='border-collapse:collapse;text-align:center;'>"
    html += "<tr>" + ''.join(f"<th>{h}</th>" for h in headers) + "</tr>"
    # Mostrar filas 0..n-1. Para cada orden j, los valores se muestran a partir de la fila j (shift hacia abajo)
    for fila in range(n):
        html += "<tr>"
        # i, xᵢ, f(xᵢ) from tabla_vis
        ix = tabla_vis[fila][0] if tabla_vis[fila][0] is not None else ''
        xi = tabla_vis[fila][1] if tabla_vis[fila][1] is not None else ''
        yi = tabla_vis[fila][2] if tabla_vis[fila][2] is not None else ''
        html += f"<td>{ix}</td><td>{xi}</td><td>{yi}</td>"
        # Ordenes: para orden j, mostrar dd[j][fila-j] si fila>=j
        for orden in range(1, max_orden+1):
            if orden < len(dd) and fila >= orden:
                idx = fila - orden
                if idx < len(dd[orden]):
                    val = dd[orden][idx]
                    html += f"<td>{'' if val is None else f'{val:.6g}'}</td>"
                else:
                    html += "<td></td>"
            else:
                html += "<td></td>"
        html += "</tr>"
    html += "</table>"
    return html

def newton_diferencias_divididas(xs, ys, x_eval):
    n = len(xs)
    if len(set(xs)) != len(xs):
        raise Exception("Hay valores de x repetidos. No se puede calcular diferencias divididas.")
    tabla, dd = tabla_diferencias_divididas(xs, ys)
    pasos = [tabla_html((tabla, dd))]
    # Coeficientes del polinomio de Newton (sin shift, diagonal principal de dd)
    coef = [dd[j][0] for j in range(n)]
    if any(c is None for c in coef):
        raise Exception("Se encontró una división por cero o datos inválidos en la tabla de diferencias divididas. Revisa que no haya x repetidos y que los datos sean correctos.")
    # Evaluación del polinomio
    resultado = coef[0]
    mult = 1
    pasos_eval = [f"P(x) = {coef[0]:.6g}"]
    for i in range(1, n):
        mult *= (x_eval - xs[i-1])
        resultado += coef[i]*mult
        pasos_eval.append(f"+ {coef[i]:.6g} * " + " * ".join([f"({x_eval} - {xs[k]})" for k in range(i)]) + f" = {coef[i]*mult:.6g}")
    pasos.extend(pasos_eval)
    pasos.append(f"Suma total: {resultado}")
    return resultado, pasos

def newton_polinomio(xs, ys):
    x = symbols('x')
    n = len(xs)
    if len(set(xs)) != len(xs):
        raise Exception("Hay valores de x repetidos. No se puede calcular diferencias divididas.")
    tabla, dd = tabla_diferencias_divididas(xs, ys)
    pasos = [tabla_html((tabla, dd))]
    coef = [dd[j][0] for j in range(n)]
    if any(c is None for c in coef):
        raise Exception("Se encontró una división por cero o datos inválidos en la tabla de diferencias divididas. Revisa que no haya x repetidos y que los datos sean correctos.")
    polinomio = coef[0]
    mult = 1
    for i in range(1, n):
        mult *= (x - xs[i-1])
        polinomio += coef[i]*mult
    pasos.append(f"Polinomio expandido: {simplify(expand(polinomio))}")
    return simplify(expand(polinomio)), pasos

@app.route('/', methods=['GET', 'POST'])
def index():
    resultado = None
    polinomio = None
    pasos = []
    pasos_poli = []
    puntos = []
    x_eval = ''
    metodo = 'lagrange'
    if request.method == 'POST':
        try:
            xs_raw = request.form.getlist('x')
            ys_raw = request.form.getlist('y')
            puntos_validos = []
            for x, y in zip(xs_raw, ys_raw):
                if x.strip() != '' and y.strip() != '':
                    try:
                        x_val = float(x)
                        y_val = float(y)
                        puntos_validos.append((x_val, y_val))
                    except Exception:
                        continue
            # Filtrar x únicos (el primero que aparece)
            xs_seen = set()
            puntos_filtrados = []
            for x, y in puntos_validos:
                if x not in xs_seen:
                    puntos_filtrados.append((x, y))
                    xs_seen.add(x)
            xs = [x for x, y in puntos_filtrados]
            ys = [y for x, y in puntos_filtrados]
            x_eval = float(request.form['x_eval'])
            metodo = request.form.get('metodo', 'lagrange')
            if metodo == 'lagrange':
                # Para Lagrange sólo mostrar la respuesta (sin pasos detallados)
                resultado, _ = lagrange_interpolacion(x_eval, xs, ys)
                polinomio, _ = lagrange_polinomio(xs, ys)
                pasos = []
                pasos_poli = []
            elif metodo == 'newton':
                resultado, pasos = newton_diferencias_divididas(xs, ys, x_eval)
                polinomio, pasos_poli = newton_polinomio(xs, ys)
            puntos = list(zip(xs, ys))
        except Exception as e:
            import traceback
            print("\n[ERROR EN FLASK]:", e)
            traceback.print_exc()
            resultado = f'Error en los datos ingresados: {e}'
    return render_template('lagrange_form.html', resultado=resultado, polinomio=polinomio, puntos=puntos, x_eval=x_eval, metodo=metodo, pasos=pasos, pasos_poli=pasos_poli)

if __name__ == '__main__':
    app.run(debug=True)
