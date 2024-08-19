### IIQ3782 - Termodinámica Avanzada


El objetivo de este curso es presentar a los estudiantes del departamento ingeniería química y bioprocesos, un conjunto interconectado de métodos computacionales necesarios en para la resolución de problemas de termodinámica avanzada. En particular, métodos que ayudan a los estudiantes a resolver problemas en las áreas equilibrios de fase utilizando ecuaciones de estado y modelos de coeficientes de actividad. 

El curso utilizará Jupyter notebooks para las experiencias (numerados del 00 al 06) para que los estudiantes practiquen sus habilidades en termodinámica computacional, usando la resolución de problemas (Ver `notebooks/`). Los notebooks están preparadas como una combinación de métodos computacionales y programación informática (en lenguaje Python) destinado a ayudar a los estudiantes para la realización de tareas y el proyecto del curso.


Los comentarios y la colaboración para mejorar este curso son bienvenidos en GitHub con `pull requests` o `issues` por correo electrónico.

Este curso utiliza Jupyter Notebooks en el lenguaje de programación Python. Se puede acceder al contenido en
las siguientes maneras:

+ La versión HTML estática de los notebooks se mostrará en el navegador actual si el notebook figura en el repositorio de código, en `notebook/`. Esto no permitirá representar siempre las fórmulas matemáticas o interactuar con el código. Alternativamente, puede renderizar los cuadernos en [NBViewer](http://nbviewer.jupyter.org/).
+ Utilice el botón verde Código arriba en la parte superior derecha de la página y descargue el repositorio a su máquina local. Descomprime el archivo. Luego use su propio servidor Jupyter Notebook (Consulte los pasos a seguir en el notebook de introducción) para navegar hasta el directorio creado por la operación de descompresión y cargar los archivos del notebook. Alternativamente, puede descargar GitHub Desktop para mantener el repositorio actualizada en su máquina local. 

Para descargar todos los archivos necesarios para completar los tutoriales, ingrese a su carpeta de trabajo y abra un terminal de comandos. Luego ingrese el siguiente comando:
```
git clone https://github.com/nfgajardo/IIQ3782.git
```

Para facilitar la instalación de los paquetes necesarios les recomendamos instalar [miniconda](https://docs.anaconda.com/free/miniconda/
). Luego, para instalar el ``conda environment`` que les instalará la mayoria de los paquetes necesarios para resolver los ``noteboks/`` tiene que hacer lo siguiente en la terminal de conda:

```
 conda env create --name iiq3782 --file=environment.yml
```

para luego activarlo usando, 

```
conda activate iiq3782
```

El resto de los paquetes tales como ``sgtpy`` y ``epcsaftpy`` se puedes instalar entrando a la carpeta correspondiente y utilizando:

```
pip install .
```

Les recomendamos mirar la documentación de cada paquete en su carpeta en caso de tener dudas con algun comando en específico. La instalación de ``phasepy`` necesita adicionalmente instalar [Visual build tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/
)

A continuación les dejamos algunos enlaces que les podrían ser utiles para comenzar con jupyter notebooks, junto con algunos paquetes de ingeniería química y termodinámicos interesantes. 

+ [Python for Chemical Engineers](https://github.com/CAChemE/Python-Chemical-Engineers)
+ [PyTherm](https://iurisegtovich.github.io/PyTherm-applied-thermodynamics/)
+ [thermo](https://github.com/CalebBell/thermo)
+ [Clapyeron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)
+ [Cheminformatics](https://github.com/PatWalters/practical_cheminformatics_tutorials)
+ [Deep learning for molecules and materials](https://dmol.pub/)
+ [Python Computations in Science and Engineering](https://kitchingroup.cheme.cmu.edu/pycse/intro.html)
+ [Machine Learning in Chemical Engineering](https://edgarsmdn.github.io/MLCE_book/intro.html)
+ [Chemical and Process Engineering Interactive Simulations](https://github.com/CAChemE/learn)
+ [Teaching and Learning with Jupyter](https://jupyter4edu.github.io/jupyter-edu-book/)
+ [Scientific Computing for Chemists with Python](https://github.com/weisscharlesj/SciCompforChemists)

Algunos cursos y artículos científicos interesantes para el proyecto
+ [Impact of Jupyter Notebook as a tool to enhance the learning process in chemical engineering modules](https://github.com/jorge-ramirez-upm/PQ-Jupyter/)
+ [Introduction to Chemical Engineering Analysis](https://github.com/jckantor/CBE20255)
+ [Computational Thermodynamics](https://kyleniemeyer.github.io/computational-thermo/content/intro.html)

Gracias de antemano por los aportes para mejorar este curso.\
Saludos,\
Equipo docente IIQ3782





