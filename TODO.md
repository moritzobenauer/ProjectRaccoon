1. Raccoon.py sortieren: Das hat alles nichts in der main verloren. / modularer aufbau
    1.1  welcome auslagern
    1.2. alle argparse gebündelt in ui

    mögliche Project struktur

raccon:
    main.py
    -src:
        -ui:
            ui kram: argpase welcome etc...
        -functions:
            __init__.py
            advanced.py
        - data:
            __init__.py
            ggf: datenstruktur für monomer schreiben?

2. globale variablen vernünftig definieren in main function

3. try excep minimieren: weshalb steht überall
 ```Python
try:
    ...
except:
    pass # raise specific error
```

4. types klären und vernünfitge docstrings schreiben


