
# IC-Eletrofisiologia-Cardiaca(português)

A simulação computacional da atividade elétrica do coração é um tema de grande interesse científico. Essas simulações buscam fornecer um melhor entendimento dos fenômenos biofísicos envolvidos e podem auxiliar no desenvolvimento de terapias para pacientes. O problema pode ser descrito por equações diferenciais parciais do tipo reação-difusão. Uma alternativa ainda pouco adotada para a solução destas equações é o uso do método de lattice Boltzmann (MLB), que tem sido cada vez mais utilizado para simular problemas da dinâmica dos fluidos. O MLB possui boas características para uso em ambientes paralelos, o que diminui o tempo de execução das simulações. Vários trabalhos implementam o método para uso em unidades de processamento gráfico, que consegue realizar várias operações aritméticas em paralelo. Porém, as implementações neste ambiente são muito dependentes do hardware utilizado. Para contornar este problema alguns trabalhos têm testado o uso do recurso chamado WebGL, onde as operações de computação gráfica exibidas no navegador são usadas para realizar computações na placa gráfica disponível. O objetivo do trabalho é a implementação do MLB utilizando WebGL para a simulação de um caso de arritmia cardíaca. Esta implementação paralela se mostrou mais eficiente do que uma implementação sequencial desenvolvida na linguagem C. Portanto, o trabalho demonstra o potencial do MLB combinado ao WebGL para a realização deste tipo de simulação.

# IC-cardiac-electrophysiology (english)

The computational simulation of the heart's electrical activity is a topic of great scientific interest. These simulations aim to provide a better understanding of the biophysical phenomena involved and can assist in the development of therapies for patients. The problem can be described by partial differential equations of the reaction-diffusion type. A still underutilized alternative for solving these equations is the use of the lattice Boltzmann method (LBM), which has been increasingly employed to simulate fluid dynamics problems. The LBM has good characteristics for use in parallel environments, which reduces the execution time of simulations. Several studies have implemented the method for use on graphics processing units (GPUs), which can perform multiple arithmetic operations in parallel. However, implementations in this environment are highly dependent on the hardware used. To address this issue, some studies have tested the use of a resource called WebGL, where the graphics computation operations displayed in the browser are used to perform computations on the available graphics card. The aim of this work is to implement the LBM using WebGL for the simulation of a case of cardiac arrhythmia. This parallel implementation has proven to be more efficient than a sequential implementation developed in the C language. Therefore, this work demonstrates the potential of combining LBM with WebGL for conducting this type of simulation.

<p align="center">
<img src="/Logos/2.png">
</p>
