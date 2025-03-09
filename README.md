# TAREA NO.3: FLUJO DE CARGA CUASI-DINAMICO CON  METODO DE PUNTO FIJO.

Este proyecto contiene una implementación en **Julia** de un flujo de potencia cuasi dinámico utilizando el método de punto fijo para resolver las tensiones nodales en un sistema eléctrico. El código incluye el cálculo estático del flujo de potencia, el procesamiento de datos de generación solar y demandas, y la simulación cuasi dinámica basada en casos típicos obtenidos mediante clustering.

---
### Estructura.

Para el desarrollo de este proyecto se dieron una serie de pasos a seguir para la implementación del flujo de carga cuasi-dinámico con método de punto fijo, los cuales se describen a continuación:

1) Calcular la Matriz de Admitancia Nodal Separada ($Y_{nn}$ y $Y_{ns}$) $\\$
Para comenzar se realizó el cálculo de la matriz de admitancia nodal o matriz $Y_{bus}$ utilizando parte de los códigos realizados en las actividades anteriores, sin embargo para el objetivo de esta actividad se separó además esta matriz en 2 submatrices $Y_{nn}$ y $Y_{ns}$, donde $Y_{nn}$ es la submatriz de admitancia entre nodos y $Y_{ns}$ es la submatriz de admitancia entre nodos y subestaciones.$\\$
Cabe resaltar que ambas matrices omiten el nodo slack, además fue necesario realizar un ajuste en la base de datos para obtener los nodos correspondientes ya que la información contenía caracteres y tipo *Strings* lo que dificultaba la indexación para la creación de las submatrices.$$\\$$

2) Calculo de Flujo Estático por Punto Fijo.$\\$
Una vez obtenida las submatrices $Y_{nn}$ y $Y_{ns}$, se procedió a realizar el cálculo del flujo de potencia estático por punto fijo basados en la expresión explicada en el marco teórico.
$$T(V_n) = Y_{nn}^{-1}\left(\left(\frac{S_n}{V_n}\right)^* - Y_{ns}V_s\right)$$
Por lo cual fué necesario la creación de una función llamada *Flujo_punto_fijo()*, en la cual se calculó previamente el vector de potencias inyectadas $S_n$ utilizando los datos aportados de la base de datos *nodes.csv*, y se inicializó un vector de tensiones nodales $V_n$ con los valores más probables en cada nodo ($V_n = 1.0 \angle0^\circ$), y se definio también la tensión en el nodo slack $V_s$ con un valor de $1.0 \angle0^\circ$; finalmente se ejecutó el flujo de carga estatico y se obtuvieron los valores de las tensiones nodales.$$\\$$

3) Procesamiento de datos de generación solar y demandas.$\\$
Las bases de datos proporcionadas para la ejecución del *flujo cuasi-dinámico* incluyen información de generación solar y demandas por intervalos de tiempo, sin embargo no existe una relación clara entre los tiempos en que se muestrean ambas bases de datos, esto es debido a que la generación solar presenta valores de geneneración en intervalos de 5 min durante todo un año, en cambio la demanda presenta valores de demanda en intervalos de 1 min durante 1 día.$\\$
Debido a esto fue necesario realizar previamente el ajuste de estos datos, comenzando por la base de datos de generación solar; para esto se convirtió esta base de datos en una matriz en la cual la cual cada columna correspondiera a un día del año, y cada fila el valor de la generación en intervalos de 5 min, de esta manera se obtuvo una matriz de 365 columnas y 288 filas.$\\$
Posteriormente se utilizó la función "K-means" para agrupar los datos de generación solar en 4 casos típicos, es decir originalmente cada día representaba un caso, sin embargo muchos de estos días presentaban valores de generación muy similares, por lo que se agruparon en 4 casos típicos.$\\$
Finalmente para la base de datos de las demandas se promedio los valores de demanda en intervalos de 5 min, de esta manera se obtuvo una base de datos con 288 filas y 55 columnas donde cada columna representaba un nodo *PQ* de los especificados en la base de datos *nodes.csv*.$$\\$$

4) Simulación cuasi-dinámica basada en casos típicos obtenidos.$\\$
Teniendo las bases de datos de demanda y generación solar ajustadas, se procedió a realizar la simulación cuasi-dinámica, para esto se creó la función *Flujo_cuasi_din()*, ya que en este caso al tener varios intervalos de tiempo en generación y demanda, el vector de potencias inyectadas $S_n$ se debía recalcular para cada intervalo de tiempo, además considerando que se tenían 4 casos típicos de generación solar, se debía recalcular el flujo de carga para cada uno de estos casos, y guardarlos de forma indepeniente para un análisis posterior; por lo cual se utilizó un *ciclo for* para recorrer cada uno de los casos típicos y otro *for anidado* para recorrer cada intervalo de tiempo, de esta forma se calcularon las tensiones nodales y se agruparon en una matriz de resultados para cada caso típico, donde cada columna representaba un intervalo de tiempo y cada fila una tensión nodal.$$\\$$

---
### Marco Teórico.

El flujo de carga es una herramienta fundamental en el análisis de sistemas eléctricos de potencia. Su propósito es determinar el estado operativo del sistema, incluyendo las tensiones nodales (*Magnitud y Ángulo*) y los flujos de potencia de las líneas, a partir de condicones conocidad como la topología del sistema, las potencias generadas y demandadas, y las características de las líneas. Este problema es no líneasl debido  a la relación entre las variables de tensión y las ecuaciones de potencia.

Las principales variables nodales son:

* $P_g$ : Potencia activa generada.

* $P_d$ : Potencia activa demandada.

* $Q_g$ : Potencia reactiva generada.

* $Q_d$ : Potencia reactiva demandada.

* $V_n$ : Tensión en el nodo n.

* $\theta_n$ : Ángulo de fase en el nodo n.

#### 1) Tipos de Nodos.

* Nodo Slack: Actúa como la referencia del sistema, con una tensión fija tipicamente de $V_s = 1.0 \angle0^\circ$. Representa la subestación principal.

* Nodos PV: Nodos generación que controlan la magnitud de la tensión y la potencia activa ($|V| y P_g$ conocidos).

* Nodos PQ: Nodos de carga o generación con potencia activa o reactiva especificadas ($P$ y $Q$ conocidos).

#### 2) Método de punto Fijo.

El método de punto fijo es una técnica iterativa para resolver ecuaciones no lineales de la forma $x=T(x)$, donde $T$ es una función contractiva. En el flujo de potencia, este método se emplea para determinar las tensiones nodales que satisfacen las ecuaciones de potencia.

##### A) Definición Matemática.

Dada una función $T: B \to B$ en un espacio de Banach $B$ (como $\mathbb{C}^2$), se dice que $T$ es una contracción si:

$$||T(x)-T(y)|| \leq k ||x-y|| \ \ con \ \ k<1$$

Por el teorema del punto fijo  de Banach, $T$ tiene un único punto fijo  $x^*$ tal que $x^* = T(x^*)$. Este punto se aproxima iterativamente mediante:

$$X^{k+1}=T(x^{k})$$

partiendo de un valor inicial $x^0$, y la secuencia {$x^{(k)}$} converge a $x^*$.

##### B) Aplicación al Flujo de Potencia.

En el flujo de potencia, el sistema se divide en:

* Nodo Slack *(s)*: Con tensión conocida $V_s$.

* Nodos No Slack *(n)*: Cuyas tensiones $V_n$ se calculan.

La matriz de admitancia nodal $[Y]$ se particiona en submatrices:
$$
\begin {bmatrix}
I_s  \\
I_n
\end{bmatrix} =
\begin {bmatrix}
Y_{ss} & Y_{sn} \\
Y_{ns} & Y_{nn}
\end{bmatrix}*
\begin {bmatrix}
V_s  \\
V_n
\end{bmatrix}
$$

Para los nodos no slack, la potencia aparente se relaciona con las corrientes mediante:

$$S_n = V_n \circ I_n^*$$

Donde  $\circ$ indica el producto de *Hadamard* (elemento a elemento).

Despejando $I_n$ y sustituyendo en la ecuación de corrientes, se obtiene:

$$ V_n = Y_{nn}^{-1}\left(\left(\frac{S_n}{V_n}\right)^* - Y_{ns}V_s\right)$$

Esta ecuación se reescribe como $V_n = T(V_n)$, donde:

$$T(V_n) = Y_{nn}^{-1}\left(\left(\frac{S_n}{V_n}\right)^* - Y_{ns}V_s\right)$$

El método de punto fijo se aplica iterativamente:

$$V_n^{(k+1)} = T(V_n^{(k)}$$

usando una estimación inicial como $V_n^{(0)} = 1 \angle0^\circ$ para todos los nodos no slack.

##### c) Criterio de Convergencia.

La convergencia se evalúa con la normal del error entre iteraciones:

$$error = ||V_n^{(k+1)}-V_n^{(k)}||$$

donde $||.||$ es la norma euclidiana. El proceso puede detenerse cuando el error cae por debajo de una tolerancia o tras un número fijo de iteraciones.

#### 3) Flujo de Potencia Cuasi Dinámico.

El flujo de potencia cuasi dinámico extiende el análisis estático para modelar variaciones temporales en las condiciones del sistema, como cambios en la generación y la demanda a lo largo del tiempo.

Para cada intervalo de tiempo se actualizan las potencias inyectadas $S_n$ según la generación y demanda del intervalo, y se calcula el flujo de potencia estatico con el método de punto fijo para obtener las tensiones nodales $V_n$.

Este enfoque permite estudiar el comportamiento del sistema bajo diferentes escenarios temporales de manera eficiente, sin recurrir a simulaciones dinámicas completas.

---

### Funciones

* **Librerias necesarias**
    - using LinearAlgebra
    - using DataFrames
    - using PrettyTables
    - using Clustering
    - using Statistics
    - using Plots
    - using CSV

* **calcular_Ynn_Yns()**
*Requiere*
    - Entradas:   
        - lines: DataFrame con los datos de las lineas del sistema
        - nodes : DataFrame con los datos de los nodos del sistema
    - Salida :    
        - Ynn : Matriz de admitancia de los nodos diferentes al Slack.
        - Yns : Matriz de admitancia de los nodos diferentes al Slack con respecto al Slack.

* **Flujo_punto_fijo()**
*Requiere*
    - Entradas: 
        - lines: DataFrame con los datos de las lineas del sistema
        - nodes : DataFrame con los datos de los nodos del sistema
    - Salida :    
        - Vn : Vector de tensiones nodales obtenido del punto fijo..
    
* **Flujo_cuasi_din()**
*Requiere*
    - Entradas: 
        - lines: DataFrame con los datos de las lineas del sistema
        - nodes : DataFrame con los datos de los nodos del sistema
    - Salida :   
        - Ninguna. Nota: La función crea K (Clusters) archivos csv con las tensiones nodales para cada día típico.
    
**Licencia**

Programa realizado por: Jean Pool Marín
jeanpool.marin@utp.edu.co
