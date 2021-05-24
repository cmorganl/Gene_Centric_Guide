import unittest

import evaluate_clusters as ec


class Tester(unittest.TestCase):
    def setUp(self) -> None:
        self.mock_ce = ec.ClusterExperiment("./")
        self.mock_ce.cluster_assignments = {"192952.BAM91895.1": ['131', '131', '132', '133'],
                                            "926571.AKB58323.1": ['12', '12'],
                                            "926571.ABS56597.1": ['1', '2']}
        return

    def tearDown(self) -> None:
        return

    def test_retrieve_lineages(self):
        acc_taxon_map = ec.retrieve_lineages([self.mock_ce])
        self.assertEqual(2, len(acc_taxon_map))
        t1 = acc_taxon_map['192952']
        t2 = acc_taxon_map['926571']
        self.assertEqual("Archaea", t1.lca(t1, t2).name)
        return

    def test_report_query_cohesiveness(self):
        completeness = self.mock_ce.report_query_completeness(['1', '2', '1', '2'])
        self.assertEqual(0, completeness)
        completeness = self.mock_ce.report_query_completeness([2, 2, 2, 2, 2])
        self.assertEqual(1.0, completeness)
        return

    def test_find_clustering_accuracy(self):
        accuracy = self.mock_ce.find_clustering_accuracy(['241', '241', '240', '238', '241', '222'])
        self.assertEqual(0.0, accuracy)
        accuracy = self.mock_ce.find_clustering_accuracy([2, 2, 2, 2, 2, 2, 2])
        self.assertEqual(1.0, accuracy)
        return

    def test_acc_line(self):
        import pandas as pd
        from plotly import express as px
        from plotly import graph_objs as go
        mock_df = pd.read_csv("acc_pd_df.csv")
        acc_plt = ec.acc_line(mock_df, palette=px.colors.qualitative.T10)
        self.assertIsInstance(acc_plt, go.Figure)
        acc_plt.show()
        return

    def test_len_bars(self):
        import pandas as pd
        from plotly import graph_objs as go
        mock_df = pd.read_csv("acc_pd_df.csv")
        bar_plt = ec.len_bars(mock_df, x_lims=[40, 101])
        self.assertIsInstance(bar_plt, go.Figure)
        bar_plt.show()
        return

    def test_smooth_sav_golay(self):
        from plotly import express as px
        df = px.data.iris()
        ret_df = ec.smooth_sav_golay(df,
                                     group_vars=["species", "species_id"],
                                     num_var="sepal_width",
                                     sort_var="sepal_length")
        self.assertEqual(150, len(ret_df))
        return

    def test_acc_summary_stats(self):
        import pandas as pd
        mock_df = pd.read_csv("acc_pd_df.csv")
        acc_summary_df = ec.acc_summary_stats(mock_df, kwargs={})
        self.assertEqual(9, len(acc_summary_df))
        return


if __name__ == '__main__':
    unittest.main()
