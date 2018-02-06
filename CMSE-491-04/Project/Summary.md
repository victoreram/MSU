### Final Project: Applying Persistent Homology to Machine Learning with Persistence Images

The following is the intro from my writeup:

  With the ongoing "data revolution", the demand for efficient methods to analyze the massive
amounts of large datasets has grown. Topology, a field that has remained relatively dormant in
providing wide-scale applications until recently, has been infused with computational techniques in
an attempt to better represent the properties of space. Computational topology has shown to be a
promising approach to analyzing noisy, complex data sets and has been used in computer vision,
medical imaging, and biological aggregation models. One of the tools used in computational topol-
ogy is persistent homology, which can describe features of interest in this dataset by representing
the homological information through a collection of points called Persistence Diagrams (PDs).

  Interest in computational topology has coincided with machine learning (ML), a prominent data
analysis technique. ML recognizes patterns among data in order to classify clusters within data and
make predictions. Techniques in ML require metrics to draw from in order to make predictions.
There's a growing interest to combine both persistence homology and machine learning in order to
refine the workflow in data analysis. It's possible to perform ML techniques using PDs and their
features as a metric when PDs are mapped into spaces which allow for it to be "fed" into a ML
algorithm.

  This project will attempt to demonstrate such a ML pipeline. Section 2 discusses the method-
ology, starting from data extraction, to data transformation, to classifying using a ML algorithm.
In this case, I used scikitlearn's linear SVM classiffier. Section 3 shows how well this model classified 
our data and is compared to results from conventional methods. Section 4 concludes with a
discussion of the limitations of this model and suggets future improvements. Much of the basis of
this project was inspired by the work of Adcock et al [PI].
