from rich.console import Console

raccoon = """  

                                                                 
                . @@@@@                   @@@@@                  
                . @@@@@@@@             @@@@@@@@                  
                  @@@@:=@@@@   @@@   @@@@##=@@@                  
                  @@@@@@::@@@@@@@@@@@@@##@+=@@@                  
                  @@@@@@########-::::::::@@@@@@                  
                .  @@@##########-::::::::::@@@+                  
        .......   @@@###     .##-::      :::@@@  ........        
   .@@@@@@@@@@   @@@#     %+    =    +*     :@@@   @@@@@@@@@@.   
   @@@@@@       @@@    @@@@@@@     @@@@@@@    @@@=      +@@@@@@  
  @@@@......   @@@  @@@@@   @@@   @@@###@@@@@  @@@@  ......@@@@  
  @@@ .....  =@@@@@@@@@@@@ .@@@===#@@%#@@@@@@@@@@@@@  ......@@@@ 
  @@@   ....   @@@@@@@@@@@@@@  ===  @@@@@@@@@@@@@@   ..   ::@@@@ 
  @@@      ...   @@@@@@+       @@@      .@@@@@@@   ....    :@@@@ 
  @@@       ....   +@@@@@.    @@@@@    .@@@@@@   ....@@@@@@@@@@@ 
  @@@        .....   -@@@@@@:       .@@@@@@@   .....@@@.  @@@@@@ 
  @@@         ...   @@@@@@@@@@@@.#@@@@@@@@@@@   ...@@@-    :@@@@ 
  @@@         ..   @@@###@@@@@@@@@@@@@@@   @@@   ...@@@@.@@@@@@@ 
  @@@         .   @@@@#####@@@@@@@@@@@.    :@@@-  ...@@@@@@@@@@@ 
  @@@         .  @@@@@@@####%@@@@@@@*    :::@@@@   ....    .@@@@ 
  @@@        ...@@@@@@@@@@####@@@@@    ::::@@@@@@   ...   ..@@@@ 
  @@@        ..@@@###@@@@@@@####@    ::::@@@...@@@  @@@@@@@@@@@@ 
  @@@        ..@@@####@@@@@@@%#    .:::@@@=....@@@   @@   ::@@@@ 
  @@@       ...@@@@@####@@@@@%#   ...@@@@....===@@   @@   .:@@@@ 
  @@@       ...@@@@@@@###@@@####-...-@@@...=====@@@  @@@@@@@@@@@ 
  @@@       ..+@@@@@@@###@@@@@%#-.:==@@@...===@@@@@  @@@@@@@@@@@ 
  @@@       ...@@@@@@@###@@@@@@@-====@@@...=@@@@@@@  .......@@@@ 
  @@@       ..+*******+++***###**++++***+++*******    . ..  @@@@ 
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  @@@                                                       @@@@ 
  @@@                         @@@@@                         @@@@ 
  @@@@                         @@@:                        =@@@  
   @@@@@                                                 @@@@@@  
   .@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
                          @@@@@@@@@@@@@                   .      
                         @@@@@@@@@@@@@@@                         
                        .@@@@@@@@@@@@@@@                         
                  .@@@@@@@@@@@@@@@@@@@@@@@@@@@.                  
                  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@                  
. ............................................................. .
      .       .....................................       .      


   
"""

creator = "Project RACCOON by Obenauer, Spauszus @ JGU Mainz 2023"
doi = "https://doi.org/10.21105/joss.06293"


def welcome():

    console = Console()
    console.print(raccoon + "\n")
    console.print(creator + "\n")
    console.print("Please cite the following references when using Project RACCOON:")
    console.print(doi + "\n", style="bold")


def tschau_kakao():
    console = Console()
    console.print("Bye Bye!" + 2 * "\n")
